import os
import sys
import json
import itertools
import sqlalchemy
import subprocess

import luigi

from luigi.contrib import sqla
from luigi.util import requires, inherits
from luigi import LocalTarget

from fieldpathogenomics.luigi.slurm import SlurmExecutableTask
from fieldpathogenomics.luigi.uv import UVExecutableTask
from fieldpathogenomics.utils import CheckTargetNonEmpty
import fieldpathogenomics.utils as utils

PIPELINE = os.path.basename(__file__).split('.')[0]
VERSION = '0.1'
luigi.auto_namespace(scope=PIPELINE)

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
        self.partition = "tgac-medium"

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
            self.partition = "tgac-medium"

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
        self.partition = "tgac-medium"

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
            self.partition = "tgac-medium"

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

# ------------------ Contiging -------------------------- #


@inherits(Trimmomatic)
class W2RapContigger(CheckTargetNonEmpty, UVExecutableTask):
    library = None
    K = luigi.IntParameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 300000
        self.n_cpu = 24
        self.host = "uv2"

    def requires(self):
        return [self.clone(Trimmomatic, library=lib.rstrip()) for lib in self.pe_libs]

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, 'contigs', 'K' + str(self.K)))

    def work_script(self):
        return '''#!/bin/bash
                    source gcc-5.2.0;
                    mkdir -p {temp_dir}
                    set -euo pipefail

                    export OMP_PROC_BIND=spread
                    export MALLOC_PER_THREAD=1
                    export w2rap='/nbi/group-data/ifs/JIC/Research-Groups/Diane-Saunders/User_Workareas/buntingd/assembly/w2rap-contigger'

                    $w2rap  --tmp_dir {temp_dir} \
                             -t {n_cpu} \
                             -m {mem} \
                             -d 16 \
                             -K {K} \
                             -o {output_dir}_temp \
                             -p pst \
                             --read_files {reads}

                  mv {output_dir}_temp {output_dir}

        '''.format(temp_dir=os.path.join(self.scratch_dir, 'pe_assembly', str(self.K)),
                   n_cpu=self.n_cpu,
                   K=self.K,
                   mem=int(0.9 * self.mem / 1000),
                   output_dir=self.output().path,
                   reads=','.join([x[0].path + ',' + x[1].path for x in self.input()]))


@requires(W2RapContigger)
class ContigStats(sqla.CopyToTable):
    columns = [
        (["K", sqlalchemy.INTEGER], {}),
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
        (["path", sqlalchemy.String(10)], {})
    ]

    connection_string = "mysql+pymysql://tgac:tgac_bioinf@tgac-db1.hpccluster/buntingd_pstgenome"
    table = "W2Rap"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_abyss(self):
        r = subprocess.run("source abyss-1.9.0; abyss-fac " + os.path.join(self.input().path, 'a.lines.fasta'),
                           stdout=subprocess.PIPE, shell=True, universal_newlines=True)
        values = r.split("\n")[1].split('\t')
        return [float(x) for x in values[:-1]] + [values[-1]]

    def rows(self):
        abyss = self.get_abyss()
        self._rows = [[self.K] + abyss]
        return self._rows

    def update_id(self):
        return hash(str(self.input().path))

# ------------------ LMP Specific -------------------------- #


@requires(FetchFastqGZ)
class UnzipFastqGZ(CheckTargetNonEmpty, SlurmExecutableTask):

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            # Set the SLURM request params for this task
            self.mem = 3000
            self.n_cpu = 1
            self.partition = "tgac-medium"

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
        self.partition = "tgac-medium"

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

                mv {cwd}/nextclip {cwd}/nextclip_done
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
        self.partition = "tgac-short"

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


@inherits(W2RapContigger)
@inherits(CompressLMP)
class MapContigs(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 750
        self.n_cpu = 8
        self.partition = "tgac-medium"

    def requires(self):
        return {'contigs': self.clone(W2RapContigger, K=200),
                'lmp': self.clone(CompressLMP)}

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, PIPELINE, VERSION, self.library, self.library + ".sam"))

    def work_script(self):
        return '''#!/bin/bash
                    source bwa-0.7.13
                    set -euo pipefail

                    #bwa index {pe_assembly}

                    bwa mem -SP -t {n_cpu} {pe_assembly} {R1} {R2} > {output}.temp

                    mv {output}.temp {output}
        '''.format(pe_assembly=os.path.join(self.input()['contigs'].path, 'a.lines.fasta'),
                   n_cpu=self.n_cpu,
                   R1=self.input()['lmp'][0].path,
                   R2=self.input()['lmp'][1].path,
                   output=self.output().path)


@requires(MapContigs)
class Sort(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 8000
        self.n_cpu = 1
        self.partition = "tgac-medium"

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


@requires(Sort)
class CollectISMetrics(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 1000
        self.n_cpu = 1
        self.partition = "tgac-short"

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

@inherits(Trimmomatic)
@inherits(CompressLMP)
class CleanedReads(luigi.WrapperTask):
    '''The LMP and the PE post-QC reads come out of different
       points in the pipeline, this task reconciles that'''

    def requires(self):
        if self.library in self.pe_libs:
            return self.clone(Trimmomatic, library=self.library)
        elif self.library in self.lmp_libs:
            return self.clone(CompressLMP, library=self.library)
        else:
            raise Exception("Unknown library")

    def output(self):
        return self.input()


@requires(CleanedReads)
class KatHist(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 4
        self.partition = "tgac-medium"

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
        self.partition = "tgac-medium"

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
        self.partition = "tgac-medium"

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
        self.mem = 4000
        self.n_cpu = 4
        self.partition = "tgac-medium"

    def requires(self):
        return {'pe': [self.clone(CleanedReads, library=lib) for lib in self.pe_libs],
                'contigs': self.clone(W2RapContigger)}

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, 'contigs', 'K' + str(self.K), self.library + "-main.mx"))

    def work_script(self):
        return '''#!/bin/bash
                    source kat-2.3.2
                    set -euo pipefail

                    kat comp -o {output_prefix} -t {n_cpu} <(zcat {reads}) {contigs}

        '''.format(output_prefix=self.output().path[:-8],
                   n_cpu=self.n_cpu,
                   reads=' '.join([x[0].path + ' ' + x[1].path for x in self.input()['pe']]),
                   contigs=os.path.join(self.input()['contigs'].path, 'a.lines.fasta'))


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
        return hists + intras + inters

# ------------------ Scaffolding -------------------------- #


@requires(W2RapContigger)
class SOAPPrep(CheckTargetNonEmpty, SlurmExecutableTask):

    prefix = luigi.Parameter(default='PST')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 28000
        self.n_cpu = 1
        self.partition = "tgac-short"

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, PIPELINE, VERSION, 'SOAP', 'K' + str(self.K), self.prefix + '.contig'))

    def work_script(self):
        return '''#!/bin/bash
                    set -euo pipefail
                    cd {cwd}
                    export soap='/usr/users/ga004/buntingd/w2rap/deps/soap_scaffolder'

                    $soap/s_prepare -g {prefix} -K 71 -c {contigs}

        '''.format(contigs=os.path.join(self.input().path, 'a.lines.fasta'),
                   cwd=self.input().path,
                   prefix=self.prefix + str(self.K))


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
        insert_sizes = {'LIB17363': 400, 'LIB26234': 800,
                        'LIB19826': 0, 'LIB19827': 0,
                        'LIB19828': 0, 'LIB19829': 0,
                        'LIB19830': 0, 'LIB19831': 0,
                        'LIB19832': 6000, 'LIB19833': 5000,
                        'LIB19834': 4000, 'LIB19835': 3500,
                        'LIB19836': 2500, 'LIB19837': 2100}

        with self.output().open('w') as fout:
            for lib in self.input()['pe'].keys():
                if insert_sizes[lib] > 0:
                    fout.write(template.format(avg_ins=insert_sizes[lib],
                                               reverse_seq='0',
                                               q1=self.input()['pe'][lib][0].path,
                                               q2=self.input()['pe'][lib][1].path))
            for lib in self.input()['lmp'].keys():
                if insert_sizes[lib] > 0:
                    fout.write(template.format(avg_ins=insert_sizes[lib],
                                               reverse_seq='1',
                                               q1=self.input()['lmp'][lib][0].path,
                                               q2=self.input()['lmp'][lib][1].path))


@inherits(SOAPPrep)
@inherits(SOAPConfig)
class SOAPMap(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 2000
        self.n_cpu = 12
        self.partition = "tgac-medium"

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, PIPELINE, VERSION, 'SOAP', 'K' + str(self.K), self.prefix + '.readOnContig.gz'))

    def requires(self):
        return {'config': self.clone(SOAPConfig),
                'contigs': self.clone(SOAPPrep)}

    def work_script(self):
        return '''#!/bin/bash
                    set -euo pipefail

                    cd {cwd}
                    export soap='/usr/users/ga004/buntingd/w2rap/deps/soap_scaffolder'

                    $soap/s_map -k 31 -s {config} -p {n_cpu} -g {prefix}

        '''.format(config=self.input()['config'].path,
                   cwd=os.path.split(self.output().path)[0],
                   n_cpu=self.n_cpu,
                   prefix=self.prefix + str(self.K))


@requires(SOAPMap)
class SOAPScaff(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 500
        self.n_cpu = 12
        self.partition = "tgac-medium"

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, PIPELINE, VERSION, 'SOAP', 'K' + str(self.K), self.prefix + '.readOnContig.gz'))

    def work_script(self):
        return '''#!/bin/bash
                    set -euo pipefail

                    cd {cwd}
                    export soap='/usr/users/ga004/buntingd/w2rap/deps/soap_scaffolder'

                    $soap/s_scaff -p {n_cpu} -g {prefix}

        '''.format(cwd=os.path.split(self.output().path)[0],
                   n_cpu=self.n_cpu,
                   prefix=self.prefix + str(self.K))


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


@inherits(PEBatchWrapper)
@inherits(LMPBatchWrapper)
@inherits(ContigStats)
@inherits(SOAPScaff)
@inherits(KATBatchWrapper)
@inherits(KatCompContigs)
class Wrapper(luigi.WrapperTask):
    library = None
    K = None

    K_list = luigi.ListParameter(default=[200])

    def requires(self):
        yield self.clone(LMPBatchWrapper)
        yield self.clone(PEBatchWrapper)
        yield self.clone(KATBatchWrapper)

        for k in self.K_list:
            yield self.clone(ContigStats, K=k)
            yield self.clone(SOAPScaff, K=k)
            for lib in self.pe_libs:
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
               '--pe-libs', json.dumps(pe_libs),
               '--lmp-libs', json.dumps(lmp_libs)] + sys.argv[3:])
