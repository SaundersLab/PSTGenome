import os
import sys
import json

import luigi

from luigi.util import requires, inherits
from luigi import LocalTarget

from fieldpathogenomics.luigi.slurm import SlurmExecutableTask
from fieldpathogenomics.luigi.uv import UVExecutableTask
from fieldpathogenomics.utils import CheckTargetNonEmpty
import fieldpathogenomics.utils as utils

PIPELINE = os.path.basename(__file__).split('.')[0]
VERSION = '0.0'
luigi.auto_namespace(scope=PIPELINE)


class FetchFastqGZ(CheckTargetNonEmpty, SlurmExecutableTask):
    '''Fetches and concatenate the fastq.gz files for ``library`` from the /reads/ server
     :param str library: library name  '''

    library = luigi.Parameter()
    base_dir = luigi.Parameter(significant=False)
    scratch_dir = luigi.Parameter(default="/tgac/scratch/buntingd/", significant=False)
    read_dir = luigi.Parameter(default="/tgac/data/reads/*DianeSaunders*", significant=False)
    lib_list = luigi.ListParameter()

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


@requires(Trimmomatic)
class KatHist(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 4
        self.partition = "tgac-medium"

    def output(self):
        return [LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, self.library, self.library + ".hist")),
                LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, self.library, self.library + ".hist.png"))]

    def work_script(self):
        return '''#!/bin/bash
                    source kat-2.3.2
                    set -euo pipefail

                    kat hist -o {output_prefix} -t {n_cpu} <(zcat {R1} {R2})

        '''.format(output_prefix=self.output()[0].path,
                   n_cpu=self.n_cpu,
                   R1=self.input()[0].path,
                   R2=self.input()[1].path)


@requires(Trimmomatic)
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


@inherits(Trimmomatic)
class KatCompInter(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 6000
        self.n_cpu = 6
        self.partition = "tgac-medium"

    def requires(self):
        return [self.clone(Trimmomatic, library=lib.rstrip()) for lib in self.lib_list]

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, self.library, self.library + ".inter-main.mx"))

    def work_script(self):
        return '''#!/bin/bash
                    source kat-2.3.2
                    set -euo pipefail

                    kat comp -n -o {output_prefix} -t {n_cpu} <(zcat {a_R1} {a_R2}) <(zcat {b_R1} {b_R2})

        '''.format(output_prefix=self.output().path[:-8],
                   n_cpu=self.n_cpu,
                   a_R1=self.input()[0][0].path,
                   a_R2=self.input()[0][1].path,
                   b_R1=self.input()[1][0].path,
                   b_R2=self.input()[1][1].path)


@inherits(Trimmomatic)
class W2RapContigger(CheckTargetNonEmpty, UVExecutableTask):
    library = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 300000
        self.n_cpu = 24
        self.host = "uv2"

    def requires(self):
        return [self.clone(Trimmomatic, library=lib.rstrip()) for lib in self.lib_list]

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, 'PE_assembly'))

    def work_script(self):
        return '''#!/bin/bash
                    source gcc-5.2.0;
                    set -euo pipefail

                    export OMP_PROC_BIND=spread
                    export MALLOC_PER_THREAD=1
                    export w2wrap='/nbi/group-data/ifs/JIC/Research-Groups/Diane-Saunders/User_Workareas/buntingd/assembly/w2rap-contigger'


                    $w2wrap  --tmp_dir {temp_dir} \
                             -t {n_cpu} \
                             -m {mem} \
                             -d 16 \
                             -o {output_dir}_temp \
                             -p run0 \
                             --read_files {reads}

                  mv {output_dir}_temp {output_dir}

        '''.format(temp_dir=os.path.join(self.scratch_dir, 'pe_assembly'),
                   n_cpu=self.n_cpu,
                   mem=int(0.9 * self.mem / 1000),
                   output_dir=self.output().path,
                   reads=','.join([x[0].path + ',' + x[1].path for x in self.input()]))


@inherits(Trimmomatic)
@inherits(W2RapContigger)
class KatCompContigs(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 4
        self.partition = "tgac-medium"

    def requires(self):
        return {'pe': self.clone(Trimmomatic),
                'contigs': self.clone(W2RapContigger)}

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, self.library, self.library + ".contigs-main.mx"))

    def work_script(self):
        return '''#!/bin/bash
                    source kat-2.3.2
                    set -euo pipefail

                    kat comp -o {output_prefix} -t {n_cpu} <(zcat {R1} {R2}) {contigs}

        '''.format(output_prefix=self.output().path[:-8],
                   n_cpu=self.n_cpu,
                   R1=self.input()['pe'][0].path,
                   R2=self.input()['pe'][1].path,
                   contigs=os.path.join(self.input()['contigs'].path, 'a.lines.fasta'))


@inherits(Trimmomatic)
@inherits(W2RapContigger)
class MapContigs(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 750
        self.n_cpu = 8
        self.partition = "tgac-medium"

    def requires(self):
        return {'contigs': self.clone(W2RapContigger),
                'reads': self.clone(Trimmomatic)}

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, self.library, self.library + ".sam"))

    def work_script(self):
        return '''#!/bin/bash
                    source bwa-0.7.13
                    set -euo pipefail

                    #bwa index {pe_assembly}

                    bwa mem -SP -t {n_cpu} {pe_assembly} {R1} {R2} > {output}.temp

                    mv {output}.temp {output}
        '''.format(pe_assembly=os.path.join(self.input()['contigs'].path, 'a.lines.fasta'),
                   n_cpu=self.n_cpu,
                   R1=self.input()['reads'][0].path,
                   R2=self.input()['reads'][1].path,
                   output=self.output().path)


@requires(MapContigs)
class Sort(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 8000
        self.n_cpu = 1
        self.partition = "tgac-short"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, self.library, self.library + ".bam"))

    def work_script(self):
        return '''#!/bin/bash
               source samtools-1.4
               set -euo pipefail

               samtools sort  --output-fmt BAM -o {output}.temp {input}

               mv {output}.temp {output}
                '''.format(input=self.input().path,
                           output=self.output().path)


@requires(Sort)
class CollectISMetrics(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 1
        self.partition = "tgac-short"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, self.library, "insert_size"))

    def work_script(self):
        return '''#!/bin/bash
               source jre-8u92
               source picardtools-2.1.1
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


@inherits(CollectISMetrics)
@inherits(KatCompInter)
@inherits(KatCompIntra)
@inherits(KatHist)
@inherits(RawFastQC)
@inherits(TrimmedFastQC)
@inherits(KatCompContigs)
@inherits(Trimmomatic)
class PerLibPipeline(luigi.WrapperTask):
    '''Wrapper task that runs all tasks on a single library'''

    def requires(self):
        return {'KatHist': self.clone(KatHist),
                'KatCompInter': self.clone(KatCompInter),
                'KatCompIntra': self.clone(KatCompIntra),
                'KatCompContigs': self.clone(KatCompContigs),
                'RawFastQC': self.clone(RawFastQC),
                'TrimmedFastQC': self.clone(TrimmedFastQC),
                'Trimmomatic': self.clone(Trimmomatic),
                'CollectISMetrics': self.clone(CollectISMetrics)}

    def output(self):
        return self.input()['Trimmomatic']


@inherits(PerLibPipeline)
class LibraryBatchWrapper(luigi.WrapperTask):
    '''Wrapper task to execute the per library part of the pipeline on all
        libraries in :param list lib_list:'''
    library = None

    def requires(self):
        return {lib: self.clone_parent(library=lib.rstrip()) for lib in self.lib_list}

    def output(self):
        return self.input()


@inherits(Trimmomatic)
class SLURMW2RapContigger(CheckTargetNonEmpty, SlurmExecutableTask):
    library = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 10000
        self.n_cpu = 12
        self.partition = "tgac-medium"

    def requires(self):
        return [self.clone(Trimmomatic, library=lib.rstrip()) for lib in self.lib_list]

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, 'PE_assembly'))

    def work_script(self):
        return '''#!/bin/bash
                    source gcc-5.2.0;
                    set -euo pipefail

                    export OMP_PROC_BIND=spread
                    export MALLOC_PER_THREAD=1
                    export w2wrap='/nbi/Research-Groups/JIC/Diane-Saunders/User_Workareas/buntingd/assembly/w2rap-contigger'


                    $w2wrap  --tmp_dir {temp_dir} \
                             -t {n_cpu} \
                             -m {mem} \
                             -d 16 \
                             -o {output_dir}_temp \
                             -p run0 \
                             --read_files {reads}

                  mv {output_dir}_temp {output_dir}

        '''.format(temp_dir=os.path.join(self.scratch_dir, 'pe_assembly'),
                   n_cpu=self.n_cpu,
                   mem=int(0.9 * self.n_cpu * self.mem / 1000),
                   output_dir=self.output().path,
                   reads=','.join([x[0].path + ',' + x[1].path for x in self.input()]))

# ----------------------------------------------------------------------- #


if __name__ == '__main__':
    os.environ['TMPDIR'] = "/tgac/scratch/buntingd"
    logger, alloc_log = utils.logging_init(log_dir=os.path.join(os.getcwd(), 'logs'),
                                           pipeline_name=PIPELINE)

    with open(sys.argv[1], 'r') as libs_file:
        lib_list = [line.rstrip() for line in libs_file if line[0] != '#']

    luigi.run(['LibraryBatchWrapper',
               '--lib-list', json.dumps(lib_list)] + sys.argv[2:])
