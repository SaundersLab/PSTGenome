import os
import sys
import json
import itertools
import sqlalchemy
import subprocess
import tempfile

import luigi

from luigi.contrib import sqla
from luigi import LocalTarget
from luigi.file import TemporaryFile

import bioluigi
from bioluigi.decorators import requires, inherits
from bioluigi.tasks import BWAIndex
from bioluigi.slurm import SlurmExecutableTask, SlurmTask
from bioluigi.utils import CheckTargetNonEmpty
import fieldpathogenomics.utils as utils

PIPELINE = os.path.basename(__file__).split('.')[0]
VERSION = '0.2'
luigi.auto_namespace(scope=PIPELINE)


# This information has to be filled in based on the output of
# CollectISMetrics, obviously hardcoding it here is not good!!!

INSERT_SIZES = {'LIB17363': 0, 'LIB26234': 800,
                'LIB19826': 0, 'LIB19827': 0,
                'LIB19828': 0, 'LIB19829': 0,
                'LIB19830': 0, 'LIB19831': 0,
                'LIB19832': 6000, 'LIB19833': 5000,
                'LIB19834': 4000, 'LIB19835': 3500,
                'LIB19836': 2500, 'LIB19837': 2100}


# ------------------ Util Tasks -------------------------- #

class AbyssFac(sqla.CopyToTable):
    columns = [
        (["Task", sqlalchemy.String(20)], {}),
        (["K", sqlalchemy.INTEGER], {}),
        (["soap_k", sqlalchemy.INTEGER], {}),
        (["n", sqlalchemy.FLOAT], {}),
        (["n:500", sqlalchemy.FLOAT], {}),
        (["L50", sqlalchemy.FLOAT], {}),
        (["min", sqlalchemy.FLOAT], {}),
        (["N80", sqlalchemy.FLOAT], {}),
        (["N50", sqlalchemy.FLOAT], {}),
        (["N20", sqlalchemy.FLOAT], {}),
        (["E-size", sqlalchemy.FLOAT], {}),
        (["max", sqlalchemy.FLOAT], {}),
        (["sum", sqlalchemy.FLOAT], {}),
        (["path", sqlalchemy.String(500)], {})
    ]

    connection_string = "mysql+pymysql://tgac:tgac_bioinf@tgac-db1.hpccluster/buntingd_pstgenome"
    table = "abyssfac"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_abyss(self):
        r = subprocess.run("source abyss-1.9.0; abyss-fac " + self.input().path,
                           stdout=subprocess.PIPE, shell=True, universal_newlines=True)
        r.check_returncode()
        values = r.stdout.split("\n")[1].split('\t')
        return [float(x) for x in values[:-1]] + [values[-1]]

    def rows(self):
        abyss = self.get_abyss()
        try:
            soap_k = self.soap_k
        except AttributeError:
            soap_k = -1
        self._rows = [[self.get_task_family()[:20]] + [self.K, soap_k] + abyss]
        return self._rows

# ------------------ Shared QC -------------------------- #


class FetchFastqGZ(CheckTargetNonEmpty, SlurmExecutableTask):
    '''Fetches and concatenate the fastq.gz files for ``library`` from the /reads/ server
     :param str library: library name  '''

    library = luigi.Parameter()
    base_dir = luigi.Parameter(significant=False)
    scratch_dir = luigi.Parameter(default="/tgac/scratch/buntingd/", significant=False)
    read_dir = luigi.Parameter(default="/tgac/data/reads/*DianeSaunders*", significant=False)

    pe_libs = luigi.ListParameter()
    lmp_libs = luigi.ListParameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 1000
        self.n_cpu = 1
        self.partition = "nbi-medium"

    def output(self):
        return [LocalTarget(os.path.join(self.scratch_dir, PIPELINE, VERSION, self.library, "raw_R1.fastq.gz")),
                LocalTarget(os.path.join(self.scratch_dir, PIPELINE, VERSION, self.library, "raw_R2.fastq.gz"))]

    def work_script(self):
        return '''#!/bin/bash -e
                  set -euo pipefail

                  find {read_dir} -name "*{library}*_R1.fastq.gz" -type f  -print | sort | xargs cat  > {R1}.temp
                  find {read_dir} -name "*{library}*_R2.fastq.gz" -type f  -print | sort | xargs cat  > {R2}.temp

                  mv {R1}.temp {R1}
                  mv {R2}.temp {R2}
                 '''.format(read_dir=self.read_dir,
                            library=self.library,
                            R1=self.output()[0].path,
                            R2=self.output()[1].path)


@requires(FetchFastqGZ)
class RawFastQC(CheckTargetNonEmpty, SlurmExecutableTask):

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            # Set the SLURM request params for this task
            self.mem = 3000
            self.n_cpu = 1
            self.partition = "nbi-medium"

        def output(self):
            return [LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, self.library, 'raw', 'R1', 'fastqc_data.txt')),
                    LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, self.library, 'raw', 'R2', 'fastqc_data.txt'))]

        def work_script(self):
            return '''#!/bin/bash
                    source fastqc-0.11.4
                    mkdir -p {output_dir}
                    set -euo pipefail

                    fastqc {R1_in} {R2_in} -o {output_dir} -t 1

                    cd {output_dir}

                    unzip raw_R1_fastqc.zip
                    sed 's/Filename\traw_R1.fastq.gz/Filename\t{lib}_R1/'  raw_R1_fastqc/fastqc_data.txt > {R1_out}.temp

                    unzip raw_R2_fastqc.zip
                    sed 's/Filename\traw_R2.fastq.gz/Filename\t{lib}_R2/'  raw_R2_fastqc/fastqc_data.txt > {R2_out}.temp

                    mv {R1_out}.temp {R1_out}
                    mv {R2_out}.temp {R2_out}
                    '''.format(output_dir=os.path.join(self.scratch_dir, PIPELINE, VERSION, self.library, 'RawFastQC'),
                               R1_in=self.input()[0].path,
                               R2_in=self.input()[1].path,
                               lib=self.library,
                               R1_out=self.output()[0].path,
                               R2_out=self.output()[1].path)

# ------------------ PE Specific -------------------------- #


@requires(FetchFastqGZ)
class Trimmomatic(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 2000
        self.n_cpu = 4
        self.partition = "nbi-medium"

    def output(self):
        return [LocalTarget(os.path.join(self.scratch_dir, PIPELINE, VERSION, self.library, "filtered_R1.fastq.gz")),
                LocalTarget(os.path.join(self.scratch_dir, PIPELINE, VERSION, self.library, "filtered_R2.fastq.gz")),
                LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, self.library, self.library + "_trimmomatic.txt"))]

    def work_script(self):
        return '''#!/bin/bash
               source jre-8u92
               source trimmomatic-0.30
               set -euo pipefail

               cd {scratch_dir}
               trimmomatic='{trimmomatic}'
               $trimmomatic PE -threads 4 {R1_in} {R2_in} -baseout temp.fastq.gz \
               ILLUMINACLIP:{adapters}:2:30:10:4 SLIDINGWINDOW:4:20 MINLEN:50 \
               2>&1 | sed 's/raw_R1.fastq.gz/{library}.fastq.gz/' > {log}.temp

               mv temp_1P.fastq.gz {R1_out}
               mv temp_2P.fastq.gz {R2_out}
               mv {log}.temp {log}

                '''.format(scratch_dir=os.path.join(self.scratch_dir, PIPELINE, VERSION, self.library),
                           trimmomatic=utils.trimmomatic.format(
                               mem=self.mem * self.n_cpu),
                           log=self.output()[2].path,
                           library=self.library,
                           R1_in=self.input()[0].path,
                           R2_in=self.input()[1].path,
                           adapters='/usr/users/ga004/buntingd/adapters.fa',
                           R1_out=self.output()[0].path,
                           R2_out=self.output()[1].path)


@requires(Trimmomatic)
class TrimmedFastQC(CheckTargetNonEmpty, SlurmExecutableTask):

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            # Set the SLURM request params for this task
            self.mem = 3000
            self.n_cpu = 1
            self.partition = "nbi-medium"

        def output(self):
            return [LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, self.library, 'trimmed', 'R1', 'fastqc_data.txt')),
                    LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, self.library, 'trimmed', 'R2', 'fastqc_data.txt'))]

        def work_script(self):
            return '''#!/bin/bash
                    source fastqc-0.11.4
                    mkdir -p {output_dir}
                    set -euo pipefail

                    fastqc {R1_in} {R2_in} -o {output_dir} -t 1

                    cd {output_dir}
                    unzip filtered_R1_fastqc.zip
                    sed 's/Filename\tfiltered_R1.fastq.gz/Filename\t{lib}_R1/'  filtered_R1_fastqc/fastqc_data.txt > {R1_out}.temp

                    unzip filtered_R2_fastqc.zip
                    sed 's/Filename\tfiltered_R2.fastq.gz/Filename\t{lib}_R2/'  filtered_R2_fastqc/fastqc_data.txt > {R2_out}.temp

                    mv {R1_out}.temp {R1_out}
                    mv {R2_out}.temp {R2_out}
                    '''.format(output_dir=os.path.join(self.scratch_dir, PIPELINE, VERSION, self.library, 'TrimmedFastQC'),
                               R1_in=self.input()[0].path,
                               R2_in=self.input()[1].path,
                               lib=self.library,
                               R1_out=self.output()[0].path,
                               R2_out=self.output()[1].path)

# ------------------ LMP Specific -------------------------- #


@requires(FetchFastqGZ)
class UnzipFastqGZ(CheckTargetNonEmpty, SlurmExecutableTask):

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            # Set the SLURM request params for this task
            self.mem = 3000
            self.n_cpu = 1
            self.partition = "nbi-short"

        def output(self):
            return [LocalTarget(os.path.join(self.scratch_dir, PIPELINE, VERSION, self.library, self.library + "_R1.fastq")),
                    LocalTarget(os.path.join(self.scratch_dir, PIPELINE, VERSION, self.library, self.library + "_R2.fastq"))]

        def work_script(self):
            return '''#!/bin/bash
                    set -euo pipefail

                    gzip -cd < {R1_in} > {R1_out}.temp
                    gzip -cd < {R2_in} > {R2_out}.temp

                    mv {R1_out}.temp {R1_out}
                    mv {R2_out}.temp {R2_out}
                    '''.format(R1_in=self.input()[0].path,
                               R2_in=self.input()[1].path,
                               R1_out=self.output()[0].path,
                               R2_out=self.output()[1].path)


@inherits(UnzipFastqGZ)
class LMP_process(CheckTargetNonEmpty, SlurmExecutableTask):

    library = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 1000
        self.n_cpu = len(self.lmp_libs)
        self.partition = "nbi-medium,RG-Diane-Saunders"

        self.cwd = os.path.join(self.scratch_dir, PIPELINE, VERSION)

    def requires(self):
        return [self.clone(UnzipFastqGZ, library=lib.rstrip()) for lib in self.lmp_libs]

    def output(self):
        return {lib: [LocalTarget(os.path.join(self.cwd, "nextclip_done", 'nextclip', lib + '_nc_ABC_R1.fastq')),
                      LocalTarget(os.path.join(self.cwd, "nextclip_done", 'nextclip', lib + '_nc_ABC_R2.fastq'))] for lib in self.lmp_libs}

    def work_script(self):
        return '''#!/bin/bash
                source flash-1.2.11;
                source python-2.7.10;
                source fastx_toolkit-0.0.13.2;
                source nextclip-1.3;

                set -euo pipefail
                export lmp_processing='/nbi/Research-Groups/JIC/Diane-Saunders/User_Workareas/buntingd/assembly/lmp_processing'

                cd {cwd}
                echo '{input}' > lmp_unzip

                echo 'Starting LMP Processing'
                echo "Using Python `which python`"

                python $lmp_processing lmp_unzip {n_cpu}

                mv -T {cwd}/nextclip {cwd}/nextclip_done
                '''.format(input='\n'.join([x[0].path + '\n' + x[1].path for x in self.input()]),
                           n_cpu=self.n_cpu,
                           cwd=self.cwd)


@requires(LMP_process)
class CompressLMP(CheckTargetNonEmpty, SlurmExecutableTask):

    library = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 1
        self.partition = "nbi-short"

    def input(self):
        return self.requires().output()[self.library]

    def output(self):
        return [LocalTarget(os.path.join(self.scratch_dir, PIPELINE, VERSION, "LMPs", self.library + '_nc_ABC_R1.fastq.gz')),
                LocalTarget(os.path.join(self.scratch_dir, PIPELINE, VERSION, "LMPs", self.library + '_nc_ABC_R2.fastq.gz'))]

    def work_script(self):
        return '''#!/bin/bash
                 set -euo pipefail

                 gzip -c < {R1_in} > {R1_out}.temp
                 gzip -c < {R2_in} > {R2_out}.temp

                 mv {R1_out}.temp {R1_out}
                 mv {R2_out}.temp {R2_out}

                 '''.format(R1_in=self.input()[0].path,
                            R2_in=self.input()[1].path,
                            R1_out=self.output()[0].path,
                            R2_out=self.output()[1].path)

# ------------------ Linked Reads Specific -------------------------- #


class LinkedReadsBAM(luigi.ExternalTask):
    library = luigi.Parameter(default='linked_reads')

    def output(self):
        return LocalTarget('/nbi/Research-Groups/JIC/Diane-Saunders/PSTGenome/longranger/pst_reads/outs/barcoded_unaligned.bam')


@requires(LinkedReadsBAM)
class BAMtoFASTQ(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 2000
        self.n_cpu = 1
        self.partition = "nbi-short"

    def output(self):
        return LocalTarget('/nbi/Research-Groups/JIC/Diane-Saunders/PSTGenome/longranger/pst_reads/outs/barcoded_unaligned.fastq.gz')

    def work_script(self):

        return '''#!/bin/bash
                    source samtools-0.1.19
                    set -euo pipefail

                    samtools view {input} |
                    perl -ne 'chomp; $line = $_;  if(/BX\:Z\:(\S{16})/){@s = split("\t", $line); print "$s[0]_$1\n$s[9]\n+\n$s[10]\n";}' |
                    gzip -c > {output}.temp

                    mv {output}.temp {output}
        '''.format(input=self.input().path,
                   output=self.output().path)

# ------------------ Cleaned Reads -------------------------- #


@inherits(Trimmomatic)
@inherits(CompressLMP)
@inherits(BAMtoFASTQ)
class CleanedReads(luigi.WrapperTask):
    '''The LMP, LR and the PE post-QC reads come out of different
       points in the pipeline, this task reconciles that'''

    def requires(self):
        if self.library in self.pe_libs:
            return self.clone(Trimmomatic, library=self.library)
        elif self.library in self.lmp_libs:
            return self.clone(CompressLMP, library=self.library)
        elif self.library == 'linked_reads':
            return self.clone(BAMtoFASTQ)
        else:
            raise Exception("Unknown library " + str(self.library))

    def output(self):
        # The Trimmomatic task is awkward and also returns the trimmomatic log
        if self.library in self.pe_libs:
            return self.input()[:2]
        else:
            return self.input()


# ------------------ Contiging -------------------------- #


@inherits(Trimmomatic)
class W2RapContigger(CheckTargetNonEmpty, SlurmExecutableTask):
    library = None
    K = luigi.IntParameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.n_cpu = 10
        self.mem = int(100000 / self.n_cpu)
        self.partition = "RG-Diane-Saunders"

    def requires(self):
        return [self.clone(Trimmomatic, library=lib.rstrip()) for lib in self.pe_libs]

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, 'contigs', 'K' + str(self.K), "a.lines.fasta"))

    def work_script(self):
        return '''#!/bin/bash
                    source gcc-5.2.0;
                    mkdir -p {temp_dir}
                    mkdir -p {output_dir}_temp
                    set -euo pipefail

                    echo `hostname`
                    echo $PATH

                    export w2rap='/nbi/group-data/ifs/JIC/Research-Groups/Diane-Saunders/User_Workareas/buntingd/assembly/w2rap-contigger'

                    $w2rap  --tmp_dir {temp_dir} \
                             -t {n_cpu} \
                             -m {mem} \
                             -d 16 \
                             -K {K} \
                             -o {output_dir}_temp \
                             -p pst \
                             --dump_all 1 \
                             --read_files {reads}

                  mv -T {output_dir}_temp {output_dir}

        '''.format(temp_dir=os.path.join(self.scratch_dir, 'pe_assembly', str(self.K)),
                   n_cpu=self.n_cpu,
                   K=self.K,
                   mem=int(0.85 * self.mem * self.n_cpu / 1000),
                   output_dir=os.path.dirname(self.output().path),
                   reads=','.join([x[0].path + ',' + x[1].path for x in self.input()]))


@inherits(CleanedReads)
class Dipspades(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 10000
        self.n_cpu = 10
        self.partition = "RG-Diane-Saunders"

    def requires(self):
        # Filter out the LMP libraries with poor insert sizes
        return {'pe': [self.clone(CleanedReads, library=lib.rstrip()) for lib in self.pe_libs],
                'lmp': [self.clone(CleanedReads, library=lib.rstrip()) for lib in self.lmp_libs if INSERT_SIZES[lib] > 0]}

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, 'contigs', 'dipspades'))

    def work_script(self):
        pe = ' '.join(["--pe{i}-1 {R1} --pe{i}-2 {R2}".format(i=i + 1, R1=x[0].path, R2=x[1].path)
                      for i, x in enumerate(self.input()['pe'])])
        lmp = ' '.join(["--mp{i}-1 {R1} --mp{i}-2 {R2}".format(i=i + 1, R1=x[0].path, R2=x[1].path)
                       for i, x in enumerate(self.input()['lmp'])])

        return '''#!/bin/bash
                  mkdir -p {output_dir}_temp
                  set -euo pipefail

                  export dipspades='/tgac/software/testing/spades/3.10.1/x86_64/bin/dipspades.py'

                  $dipspades -o {output_dir}_temp \
                             --thread {n_cpu} \
                             --memory {mem} \
                              {pe} \
                              {lmp}

                  mv -T {output_dir}_temp {output_dir}
        '''.format(n_cpu=self.n_cpu,
                   mem=int(0.9 * self.mem * self.n_cpu / 1000),
                   output_dir=self.output().path,
                   pe=pe,
                   lmp=lmp)


@requires(W2RapContigger)
class ContigStats(AbyssFac):
    pass


@requires(W2RapContigger)
class BWAIndexContigs(BWAIndex):
    pass


# ------------------ LMP Insert Sizes -------------------------- #


@inherits(BWAIndexContigs)
@inherits(CompressLMP)
class MapContigs(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 4
        self.partition = "nbi-medium"

    def requires(self):
        return {'contigs': self.clone(BWAIndexContigs, K=200),
                'lmp': self.clone(CompressLMP)}

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, PIPELINE, VERSION, self.library, self.library + ".sam"))

    def work_script(self):
        return '''#!/bin/bash
                    source bwa-0.7.13
                    set -euo pipefail

                    bwa mem -SP -t {n_cpu} {pe_assembly} {R1} {R2} > {output}.temp

                    mv {output}.temp {output}
        '''.format(pe_assembly=os.path.join(self.input()['contigs'].path[:-4]),
                   n_cpu=self.n_cpu,
                   R1=self.input()['lmp'][0].path,
                   R2=self.input()['lmp'][1].path,
                   output=self.output().path)


@requires(MapContigs)
class SortLMP(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 8000
        self.n_cpu = 1
        self.partition = "nbi-medium"

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, PIPELINE, VERSION, self.library, self.library + ".bam"))

    def work_script(self):
        return '''#!/bin/bash
               source samtools-1.4
               set -euo pipefail

               samtools sort --output-fmt BAM -o {output}.temp {input}

               mv {output}.temp {output}
                '''.format(input=self.input().path,
                           output=self.output().path)


@requires(SortLMP)
class CollectISMetrics(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 1000
        self.n_cpu = 1
        self.partition = "nbi-short"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, self.library, "insert_size"))

    def work_script(self):
        return '''#!/bin/bash
               source jre-8u92
               source picardtools-2.1.1
               source R-3.3.1;
               picard='{picard}'
               set -euo pipefail

               $picard CollectInsertSizeMetrics \
                       VERBOSITY=ERROR \
                       QUIET=true \
                       I={input} \
                       O={output}.temp \
                       M=0.5 \
                       H={output}.pdf

               mv {output}.temp {output}
                '''.format(input=self.input().path,
                           output=self.output().path,
                           picard=utils.picard.format(mem=self.mem * self.n_cpu))

# ------------------ Kmer Spectra -------------------------- #


@requires(CleanedReads)
class KatHist(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 4
        self.partition = "nbi-medium,tgac-medium,RG-Diane-Saunders"

    def output(self):
        return [LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, self.library, self.library + ".hist"))]

    def work_script(self):
        return '''#!/bin/bash
                    source kat-2.3.2
                    set -euo pipefail

                    kat hist -o {output_prefix} -t {n_cpu} <(zcat {R1} {R2})

        '''.format(output_prefix=self.output()[0].path,
                   n_cpu=self.n_cpu,
                   R1=self.input()[0].path,
                   R2=self.input()[1].path)


@requires(CleanedReads)
class KatCompIntra(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 6000
        self.n_cpu = 6
        self.partition = "nbi-medium,tgac-medium"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, self.library, self.library + ".intra-main.mx"))

    def work_script(self):
        return '''#!/bin/bash
                    source kat-2.3.2
                    set -euo pipefail

                    kat comp -n -o {output_prefix} -t {n_cpu} <(zcat {R1}) <(zcat {R2})

        '''.format(output_prefix=self.output().path[:-8],
                   n_cpu=self.n_cpu,
                   R1=self.input()[0].path,
                   R2=self.input()[1].path)


@inherits(CleanedReads)
class KatCompInter(CheckTargetNonEmpty, SlurmExecutableTask):
    libA = luigi.Parameter()
    libB = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 6000
        self.n_cpu = 6
        self.partition = "nbi-medium,tgac-medium"

    def requires(self):
        return [self.clone(CleanedReads, library=self.libA),
                self.clone(CleanedReads, library=self.libB)]

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, 'KAT', "{0}_{1}-main.mx".format(self.libA, self.libB)))

    def work_script(self):
        return '''#!/bin/bash
                    source kat-2.3.2
                    set -euo pipefail

                    kat comp -n -o {output_prefix} -t {n_cpu} <(zcat {a_R1} {a_R2}) <(zcat {b_R1} {b_R2})
                    kat plot spectra-mx -o {output_prefix}.spectra_mx.png -x 100 --intersection {output_prefix}-main.mx

        '''.format(output_prefix=self.output().path[:-8],
                   n_cpu=self.n_cpu,
                   a_R1=self.input()[0][0].path,
                   a_R2=self.input()[0][1].path,
                   b_R1=self.input()[1][0].path,
                   b_R2=self.input()[1][1].path)


@inherits(CleanedReads)
@inherits(W2RapContigger)
class KatCompContigs(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 6000
        self.n_cpu = 6
        self.partition = "nbi-medium,tgac-medium"

    def requires(self):
        return {'reads': self.clone(CleanedReads, library=self.library),
                'contigs': self.clone(W2RapContigger)}

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, 'contigs', 'K' + str(self.K), self.library + "-main.mx"))

    def work_script(self):
        reads = (self.input()['reads'][0].path + ' ' + self.input()['reads'][1].path
                 if isinstance(self.input()['reads'], list)
                 else self.input()['reads'].path)
        return '''#!/bin/bash
                    source kat-2.3.2
                    set -euo pipefail

                    kat comp -o {output_prefix} -t {n_cpu} <(zcat {reads}) {contigs}

        '''.format(output_prefix=self.output().path[:-8],
                   n_cpu=self.n_cpu,
                   reads=reads,
                   contigs=self.input()['contigs'].path)


@inherits(KatHist)
@inherits(KatCompIntra)
@inherits(KatCompInter)
class KATBatchWrapper(luigi.WrapperTask):
    '''Wrapper task to execute the per library part of the pipeline on all
        libraries in :param list lib_list:'''
    library = None
    libA = None
    libB = None

    def requires(self):
        libs = self.pe_libs + self.lmp_libs
        hists = [self.clone(KatHist, library=lib) for lib in libs]
        intras = [self.clone(KatCompIntra, library=lib) for lib in libs]
        inters = [self.clone(KatCompInter, libA=a, libB=b) for a, b in itertools.combinations(libs, 2)]
        return hists #+ intras + inters

# ------------------ Homozygous filter -------------------------- #


@requires(W2RapContigger)
class Redundans(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 4
        self.partition = "nbi-medium,RG-Diane-Saunders"

    def output(self):
        return [LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, 'contigs', 'redundans', str(self.K), "homozygous.fa")),
                LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, 'contigs', 'redundans', str(self.K), "heterozygous.fa"))]

    def work_script(self):

        return '''#!/bin/bash
                    source redundans-0.13;
                    set -euo pipefail

                    redundans -f {contigs} \
                              -o {temp_dir} \
                              --identity 0.65 \
                              --overlap 0.85  \
                              --minLength 300 \
                              --noscaffolding --nogapclosing  -v -t{n_cpu}

                    mv {temp_dir}/contigs.reduced.fa {output_hom}
                    mv {temp_dir}/contigs.reduced.fa.hetero.tsv {output_het}
        '''.format(contigs=self.input().path,
                   n_cpu=self.n_cpu,
                   temp_dir=os.path.join(tempfile.mkdtemp(), 'redundans'),
                   output_hom=self.output()[0].path,
                   output_het=self.output()[1].path)


@requires(Redundans)
class HomozygousContigs(luigi.WrapperTask):
    def output(self):
        return self.input()[0]


@requires(HomozygousContigs)
class BWAIndexContigs(BWAIndex):
    pass

# ------------------ Scaffolding -------------------------- #


@requires(HomozygousContigs)
class SOAPPrep(CheckTargetNonEmpty, SlurmExecutableTask):

    prefix = luigi.Parameter(default='PST')
    soap_k = luigi.IntParameter(default=71)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 28000
        self.n_cpu = 1
        self.partition = "nbi-short"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, 'SOAP', 'K' + str(self.K), 'k' + str(self.soap_k), self.prefix + '.contig'))

    def work_script(self):
        return '''#!/bin/bash
                    set -euo pipefail
                    cd {cwd}
                    export soap='/usr/users/ga004/buntingd/w2rap/deps/soap_scaffolder'

                    $soap/s_prepare -g {prefix} -K {k} -c {contigs}

        '''.format(contigs=self.input().path,
                   cwd=os.path.dirname(self.output().path),
                   prefix=self.output().path[:-7],
                   k=self.soap_k)


@inherits(LMP_process)
@inherits(Trimmomatic)
class SOAPConfig(CheckTargetNonEmpty, luigi.Task):

    def requires(self):
        return {'lmp': self.clone(LMP_process),
                'pe': {lib: self.clone(Trimmomatic, library=lib) for lib in self.pe_libs}}

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, 'SOAP', 'config.txt'))

    def run(self):
        template = '[LIB]\navg_ins={avg_ins}\nreverse_seq={reverse_seq}\nq1={q1}\nq2={q2}\n'

        with self.output().open('w') as fout:
            for lib in self.input()['pe'].keys():
                if INSERT_SIZES[lib] > 0:
                    fout.write(template.format(avg_ins=INSERT_SIZES[lib],
                                               reverse_seq='0',
                                               q1=self.input()['pe'][lib][0].path,
                                               q2=self.input()['pe'][lib][1].path))
            for lib in self.input()['lmp'].keys():
                if INSERT_SIZES[lib] > 0:
                    fout.write(template.format(avg_ins=INSERT_SIZES[lib],
                                               reverse_seq='1',
                                               q1=self.input()['lmp'][lib][0].path,
                                               q2=self.input()['lmp'][lib][1].path))


@inherits(SOAPPrep)
@inherits(SOAPConfig)
class SOAPMap(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 5
        self.partition = "nbi-medium,RG-Diane-Saunders,tgac-medium"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, 'SOAP', 'K' + str(self.K), 'k' + str(self.soap_k), self.prefix + '.readOnContig.gz'))

    def requires(self):
        return {'config': self.clone(SOAPConfig),
                'contigs': self.clone(SOAPPrep)}

    def work_script(self):
        return '''#!/bin/bash
                    set -euo pipefail

                    cd {cwd}
                    export soap='/usr/users/ga004/buntingd/w2rap/deps/soap_scaffolder'

                    $soap/s_map -k {k} -s {config} -p {n_cpu} -g {contigs}

        '''.format(config=self.input()['config'].path,
                   cwd=os.path.dirname(self.output().path),
                   n_cpu=self.n_cpu,
                   k=self.soap_k,
                   contigs=self.input()['contigs'].path[:-7])


@requires(SOAPMap)
class SOAPScaff(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 2000
        self.n_cpu = 4
        self.partition = "nbi-medium,RG-Diane-Saunders"

    def output(self):
        return [LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, 'SOAP', 'K' + str(self.K), 'k' + str(self.soap_k), self.prefix + '.scafSeq')),
                LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, 'SOAP', 'K' + str(self.K), 'k' + str(self.soap_k), self.prefix + '.scaf'))]

    def work_script(self):
        return '''#!/bin/bash
                    set -euo pipefail

                    cd {cwd}
                    export soap='/usr/users/ga004/buntingd/w2rap/deps/soap_scaffolder'

                    $soap/s_scaff -p {n_cpu} -g {prefix}

        '''.format(cwd=os.path.dirname(self.output()[0].path),
                   n_cpu=self.n_cpu,
                   prefix=self.input().path[:-16])


@requires(SOAPScaff)
class SOAPNremap(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 6000
        self.n_cpu = 1
        self.partition = "nbi-medium"
        self.rm_tmp = False

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, 'SOAP' + str(self.soap_k), 'K' + str(self.K), self.prefix + str(self.K) + '.scaf.fasta'))

    def work_script(self):
        soap_base = self.input()[0].path[:-8]

        return '''#!/bin/bash
                    source python-2.7.11
                    set -euo pipefail

                    export SOAP_n_remapper='/usr/users/ga004/buntingd/w2rap/SOAP_n_remapper.py'

                    /tgac/software/testing/python/2.7.11/x86_64/bin/python $SOAP_n_remapper {contig_pos_in_scaff} {scaffolds_file} {contigs_file} {output}.temp

                    mv {output}.temp {output}
        '''.format(contig_pos_in_scaff=soap_base + '.contigPosInscaff',
                   scaffolds_file=soap_base + '.scafSeq',
                   contigs_file=soap_base + '.contig',
                   output=self.output().path)


@requires(SOAPNremap)
class ScaffoldStats(AbyssFac):
    pass


@requires(SOAPNremap)
class BWAIndexScaffolds(BWAIndex):
    pass

# ------------------ Gap Filling -------------------------- #


@inherits(CleanedReads)
class AbyssBloomBuild(CheckTargetNonEmpty, SlurmExecutableTask):

    bloom_k = luigi.IntParameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 750
        self.n_cpu = 20
        self.partition = "nbi-medium"

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, PIPELINE, VERSION, "bloomfilters", str(self.bloom_k) + ".bloom"))

    def requires(self):
        return [self.clone(CleanedReads, library=lib) for lib in list(self.lmp_libs) + list(self.pe_libs) + ['linked_reads']]

    def work_script(self):
        input = [' '.join([x.path for x in y]) if isinstance(y, list) else y.path for y in self.input()]
        return '''#!/bin/bash
                    source abyss-2.0.2;
                    set -euo pipefail

                    abyss-bloom build -k{k} -b{size} -j{n_cpu} {output}.temp {input}

                    mv {output}.temp {output}
        '''.format(k=self.bloom_k,
                   size='10G',
                   input=' '.join(input),
                   n_cpu=self.n_cpu,
                   output=self.output().path)


@inherits(SOAPNremap)
@inherits(AbyssBloomBuild)
class AbyssSealer(CheckTargetNonEmpty, SlurmExecutableTask):

    sealer_klist = luigi.ListParameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 3000
        self.n_cpu = 4
        self.partition = "nbi-long"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, "sealer", 'SOAP' + str(self.soap_k), 'K' + str(self.K), 'K' + str(self.K) + "_scaffold.fa"))

    def requires(self):
        return {'bloomfilters': [self.clone(AbyssBloomBuild, bloom_k=k) for k in self.sealer_klist],
                'scaffolds': self.clone(SOAPNremap)}

    def work_script(self):
        return '''#!/bin/bash
                    source abyss-2.0.2;
                    mkdir -p {output}/temp
                    set -euo pipefail

                    abyss-sealer {k_args} \
                                 --max-paths=25 \
                                 --max-gap-length=6000 \
                                 --max-branches=10000  \
                                 --verbose \
                                 --flank-length=220 \
                                 -j {n_cpu} \
                                 -o {output}/temp/{prefix} \
                                 -S {scaffolds} {bloomfilters}

                    mv {output}/temp/{prefix}* {output}/
        '''.format(k_args=' '.join(['-k' + str(k) for k in self.sealer_klist]),
                   bloomfilters=' '.join(['-i ' + x.path for x in self.input()['bloomfilters']]),
                   scaffolds=self.input()['scaffolds'].path,
                   n_cpu=self.n_cpu,
                   output=os.path.dirname(self.output().path),
                   prefix='K' + str(self.K))


@requires(AbyssSealer)
class ScaffoldToContigs(CheckTargetNonEmpty, SlurmTask):
    '''Gap filling improves contiguity so split the gap filled scaffolds into longer contigs'''
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 1000
        self.n_cpu = 1
        self.partition = "nbi-short"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, "sealer", 'SOAP' + str(self.soap_k), 'K' + str(self.K), 'K' + str(self.K) + "_contigs.fa"))

    def work(self):
        import Bio.SeqIO
        import re

        N_min = 5  # Split on runs of N's greater than 5
        r = re.compile('N{' + str(N_min) + ',}')

        with self.input().open('r') as fin, self.output().open('w') as fout:
            scaffolds = Bio.SeqIO.parse(fin, 'fasta')

            contig_n = 0
            for scaf in scaffolds:
                for cont in r.split(str(scaf.seq)):
                    Bio.SeqIO.write(Bio.SeqRecord.SeqRecord(id='contig_{0} {1}'.format(contig_n, scaf.id), seq=Bio.Seq.Seq(cont)),
                                    fout, 'fasta')


@requires(ScaffoldToContigs)
class ContigsGFStats(AbyssFac):
    pass


@requires(AbyssSealer)
class BWAIndexGFScaffolds(BWAIndex):
    pass
# ------------------ Linked Reads  -------------------------- #


class LinkedReadsBAM(luigi.ExternalTask):

    def output(self):
        return LocalTarget('/nbi/Research-Groups/JIC/Diane-Saunders/PSTGenome/longranger/pst_reads/outs/barcoded_unaligned.bam')


@requires(LinkedReadsBAM)
class BAMtoFASTQ(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 2000
        self.n_cpu = 1
        self.partition = "nbi-medium"

    def output(self):
        return LocalTarget('/nbi/Research-Groups/JIC/Diane-Saunders/PSTGenome/longranger/pst_reads/outs/barcoded_unaligned.fastq.gz')

    def work_script(self):

        return '''#!/bin/bash
                    source samtools-0.1.19
                    set -euo pipefail

                    samtools view {input} |
                    perl -ne 'chomp; $line = $_;  if(/BX\:Z\:(\S{16})/){@s = split("\t", $line); print "$s[0]_$1\n$s[9]\n+\n$s[10]\n";}' |
                    gzip -c > {output}.temp

                    mv {output}.temp {output}
        '''.format(input=self.input().path,
                   output=self.output().path)


@requires(reads=BAMtoFASTQ, contigs=BWAIndexScaffolds)
class MapLinkedReadsScaffolds(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 500
        self.n_cpu = 24
        self.partition = "nbi-medium"

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, PIPELINE, VERSION, "linked_reads", 'SOAP' + str(self.soap_k), "K" + str(self.K), "barcoded_aligned.bam"))

    def work_script(self):
        return '''#!/bin/bash
                    source bwa-0.7.13
                    source samtools-1.4
                    set -euo pipefail

                    cd {cwd}

                    bwa mem -t {n_cpu} -p {pe_assembly} {reads} |
                    samtools view -Sb - |
                    samtools sort -n - > {output}.temp

                    mv {output}.temp {output}
        '''.format(pe_assembly=os.path.join(self.input()['contigs'].path[:-4]),
                   n_cpu=self.n_cpu - 2,
                   reads=self.input()['reads'].path,
                   output=self.output().path,
                   cwd=os.path.dirname(self.output().path))


@requires(reads=BAMtoFASTQ, contigs=BWAIndexGFScaffolds)
class MapLinkedReadsGFScaffolds(MapLinkedReadsScaffolds):
    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, PIPELINE, VERSION, "sealer", "linked_reads", 'SOAP' + str(self.soap_k), "K" + str(self.K), "barcoded_aligned.bam"))


# ------------------ Architect Scaffolding  -------------------------- #

@requires(MapLinkedReadsScaffolds)
class BamToContainment(CheckTargetNonEmpty, SlurmTask):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 8000
        self.n_cpu = 1
        self.partition = "nbi-short"

        self.threshold = 10  # The minimum number of reads per barcode

    def output(self):
        return LocalTarget(self.input().path[:-4] + ".containment")

    def work(self):
        import pysam

        contig_interval_map = {}
        contig_well_suport_counts = {}
        barcodes_to_well = {}
        n = 0

        samfile = pysam.Samfile(self.input().path, "rb")

        for read in samfile:
            if (not read.positions) or read.mapping_quality < 30:
                continue
            bc = read.qname.split('_')[1]

            # Number the barcodes and use this as the well id
            if bc not in barcodes_to_well:
                barcodes_to_well[bc] = len(barcodes_to_well.keys())
            well = barcodes_to_well[bc]

            contig = samfile.getrname(read.tid)
            start, end = min(read.positions), max(read.positions)

            if contig not in contig_interval_map:
                contig_interval_map[contig] = dict()
                contig_well_suport_counts[contig] = dict()
            if well not in contig_interval_map[contig]:
                contig_interval_map[contig][well] = [[start], [end]]
                contig_well_suport_counts[contig][well] = 0

            start_positions, end_positions = contig_interval_map[contig][well]

            if start < start_positions[-1]:
                start_positions.append(start)
                start_positions.sort()
                start_positions = start_positions[:10]
            if end > end_positions[0]:
                end_positions.append(end)
                end_positions.sort()
                end_positions = end_positions[1:]

            contig_interval_map[contig][well] = start_positions, end_positions
            contig_well_suport_counts[contig][well] += 1

            n += 1
            if n % 1000000 == 0:
                print('reads processed: {0}'.format(n))

        # write containment
        with self.output().open('w') as c:
            for contig, interval_map in contig_interval_map.items():
                for well, interval in interval_map.items():
                    if contig_well_suport_counts[contig][well] > self.threshold:
                        start, end = interval[0][-1], interval[1][0]
                        c.write('W\t%s\t%d\t%d\t%d\n' % (contig, well, start, end))


@requires(fasta=SOAPNremap, containment=BamToContainment)
class ArchitectScaffolds(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 8000
        self.n_cpu = 1
        self.partition = "nbi-medium"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, 'architect', 'SOAP' + str(self.soap_k), "K" + str(self.K), 'architect'))

    def work_script(self):

        return '''#!/bin/bash
               source /usr/users/ga004/buntingd/FP_dev/misc/python2.7/bin/activate
               set -euo pipefail

                python /usr/users/ga004/buntingd/architect/architect.py scaffold \
                --fasta {fasta} \
                --containment {containment} \
                --out {output} \
                #--edges EDGES]


                '''.format(fasta=self.input()['fasta'].path,
                           containment=self.input()['containment'].path,
                           output=self.output().path)


@requires(ArchitectScaffolds)
class ArchitectStats(AbyssFac):
    pass

# ------------------ ARCS Scaffolding  -------------------------- #


@requires(contigs=SOAPNremap, reads=MapLinkedReadsScaffolds)
class ArcsLinksScaffolds(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 8000
        self.n_cpu = 1
        self.partition = "nbi-medium"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, 'linked_reads', 'SOAP' + str(self.soap_k), "K" + str(self.K), "K" + str(self.K) + ".scaffolds.fa"))

    def work_script(self):
        self.temp = TemporaryFile()

        return '''#!/bin/bash
               source ARCS-1.0.0
               source LINKS-1.8.5
               mkdir -p {output}_temp
               set -euo pipefail

               cd {output}_temp
               echo '{reads}' > {temp}

               echo 'Starting ARCS'
               arcs -v -f {contigs} -a {temp} -b {prefix}

               echo 'Starting makeTSVfile'
               makeTSVfile.py {prefix}_original.gv {prefix}.tigpair_checkpoint.tsv {contigs}

               echo 'Starting LINKS'
               perl /tgac/software/testing/LINKS/1.8.5/x86_64/LINKS -f {contigs} \
                                                                    -s /dev/null \
                                                                    -k 20 \
                                                                    -b {prefix} \
                                                                    -l 5 \
                                                                    -t 2 \
                                                                    -a 0.7
                # The ARCS paper recommends a = 0.7-0.9
                # The arcs -e parameter is the length of the ends, should be smaller for smaller conitgs
               echo 'done'
               mv -T {output}_temp {output}
                '''.format(contigs=self.input()['contigs'].path,
                           reads=self.input()['reads'].path,
                           output=os.path.dirname(self.output().path),
                           temp=self.temp.path,
                           prefix="K" + str(self.K))


@requires(ArcsLinksScaffolds)
class ARCSStats(AbyssFac):
    pass


@requires(contigs=AbyssSealer, reads=MapLinkedReadsGFScaffolds)
class ArcsLinksGFScaffolds(ArcsLinksScaffolds):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 8000
        self.n_cpu = 1
        self.partition = "nbi-medium"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, 'sealer', 'linked_reads', 'SOAP' + str(self.soap_k), "K" + str(self.K), "K" + str(self.K) + ".scaffolds.fa"))


@requires(ArcsLinksGFScaffolds)
class ARCSGFStats(AbyssFac):
    pass
# ------------------ Supernova -------------------------- #


class SupernovaMegabubbles(luigi.ExternalTask):
    def output(self):
        return LocalTarget("/nbi/Research-Groups/JIC/Diane-Saunders/PSTGenome/supernova/1.2.0/pst88/fasta/megabubbles.fasta")


@inherits(SupernovaMegabubbles)
@inherits(CleanedReads)
class KatCompSupernova(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 6000
        self.n_cpu = 6
        self.partition = "nbi-medium,tgac-medium"

    def requires(self):
        return {'reads': self.clone(CleanedReads, library=self.library),
                'contigs': self.clone(SupernovaMegabubbles)}

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, 'supernova', 'KAT', self.library + "-main.mx"))

    def work_script(self):
        reads = (self.input()['reads'][0].path + ' ' + self.input()['reads'][1].path
                 if isinstance(self.input()['reads'], list)
                 else self.input()['reads'].path)
        return '''#!/bin/bash
                    source kat-2.3.2
                    set -euo pipefail

                    kat comp -o {output_prefix} -t {n_cpu} <(zcat {reads}) {contigs}

        '''.format(output_prefix=self.output().path[:-8],
                   n_cpu=self.n_cpu,
                   reads=reads,
                   contigs=self.input()['contigs'].path)


@requires(SupernovaMegabubbles)
class BWAIndexMegabubbles(BWAIndex):
    pass


@inherits(BWAIndexMegabubbles)
@inherits(W2RapContigger)
class MapContigsMegabubbles(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 200
        self.n_cpu = 4
        self.partition = "nbi-medium"

    def requires(self):
        return {'contigs': self.clone(BWAIndexMegabubbles),
                'reads': self.clone(W2RapContigger)}

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, "supernova", 'megabubbles', "K" + str(self.K), "mapped.bam"))

    def work_script(self):
        return '''#!/bin/bash
                    source bwa-0.7.13
                    source samtools-1.4
                    set -euo pipefail

                    cd {cwd}

                    bwa mem -t {n_cpu} -p {pe_assembly} {reads} |
                    samtools view -Sb - >{output}.temp

                    mv {output}.temp {output}
        '''.format(pe_assembly=os.path.join(self.input()['contigs'].path[:-4]),
                   n_cpu=self.n_cpu - 1,
                   reads=self.input()['reads'].path,
                   output=self.output().path,
                   cwd=os.path.dirname(self.output().path))


class SupernovaPseudoHaps(luigi.ExternalTask):
    def output(self):
        return [LocalTarget("/nbi/Research-Groups/JIC/Diane-Saunders/PSTGenome/supernova/pst88/fasta/pseudohap2.1.fasta.gz"),
                LocalTarget("/nbi/Research-Groups/JIC/Diane-Saunders/PSTGenome/supernova/pst88/fasta/pseudohap2.2.fasta.gz")]

# ------------------ Wrapper Tasks -------------------------- #


@inherits(RawFastQC)
@inherits(TrimmedFastQC)
@inherits(Trimmomatic)
class PEPerLibPipeline(luigi.WrapperTask):
    '''Wrapper task that runs all tasks on a single library'''

    def requires(self):
        return {'RawFastQC': self.clone(RawFastQC),
                'TrimmedFastQC': self.clone(TrimmedFastQC),
                'Trimmomatic': self.clone(Trimmomatic)}

    def output(self):
        return self.input()['Trimmomatic']


@inherits(PEPerLibPipeline)
class PEBatchWrapper(luigi.WrapperTask):
    '''Wrapper task to execute the per library part of the pipeline on all
        libraries in :param list lib_list:'''
    library = None

    def requires(self):
        return [self.clone(PEPerLibPipeline, library=lib.rstrip()) for lib in self.pe_libs]

    def output(self):
        return self.input()


@inherits(RawFastQC)
@inherits(CollectISMetrics)
@inherits(CompressLMP)
class LMPPerLibPipeline(luigi.WrapperTask):
    '''Wrapper task that runs all tasks on a single library'''

    def requires(self):
        return [self.clone(RawFastQC),
                self.clone(CompressLMP),
                self.clone(CollectISMetrics)]

    def output(self):
        return self.input()


@inherits(LMPPerLibPipeline)
class LMPBatchWrapper(luigi.WrapperTask):
    '''Wrapper task to execute the per library part of the pipeline on all
        libraries in :param list lib_list:'''
    library = None

    def requires(self):
        return [self.clone(LMPPerLibPipeline, library=lib.rstrip()) for lib in self.lmp_libs]

    def output(self):
        return self.input()


@inherits(PEBatchWrapper, LMPBatchWrapper, ContigStats,
          KATBatchWrapper, KatCompContigs, Dipspades,
          ARCSStats, ARCSGFStats, MapContigsMegabubbles,
          ContigsGFStats, ArchitectStats)
class Wrapper(luigi.WrapperTask):
    library = None
    K = None
    bloom_k = None
    bloom_size = None
    K_list = luigi.ListParameter(default=[200])
    soap_klist = luigi.ListParameter()

    def requires(self):
        yield self.clone(LMPBatchWrapper)
        yield self.clone(PEBatchWrapper)
        yield self.clone(KATBatchWrapper)
        yield self.clone(MapContigsMegabubbles, K=280)

        yield self.clone(Dipspades)

        for lib in list(self.pe_libs):
            yield self.clone(KatCompSupernova, library=lib)

        for k in self.K_list:
            yield self.clone(ContigStats, K=k)

            for soap_k in self.soap_klist:
                yield self.clone(ScaffoldStats, K=k, soap_k=soap_k)
                yield self.clone(ARCSStats, K=k, soap_k=soap_k)
                yield self.clone(ARCSGFStats, K=k, soap_k=soap_k)
                yield self.clone(ContigsGFStats, K=k, soap_k=soap_k)
                yield self.clone(ArchitectStats, K=k, soap_k=soap_k)

            for lib in list(self.pe_libs):
                yield self.clone(KatCompContigs, K=k, library=lib)

# ----------------------------------------------------------------------------------------------------- #


if __name__ == '__main__':
    os.environ['TMPDIR'] = "/tgac/scratch/buntingd"
    logger, alloc_log = utils.logging_init(log_dir=os.path.join(os.getcwd(), 'logs'),
                                           pipeline_name=PIPELINE)

    with open(sys.argv[1], 'r') as pe_libs_file:
        pe_libs = [line.rstrip() for line in pe_libs_file if line[0] != '#']

    with open(sys.argv[2], 'r') as lmp_libs_file:
        lmp_libs = [line.rstrip() for line in lmp_libs_file if line[0] != '#']

    luigi.run(['Wrapper',
               '--sealer-klist', json.dumps([200, 180, 160, 140]),
               '--soap-klist', json.dumps([21, 31, 51, 71, 91, ]),
               '--pe-libs', json.dumps(pe_libs),
               '--lmp-libs', json.dumps(lmp_libs)] + sys.argv[3:])
