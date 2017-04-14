import os
import sys
import json

import luigi

from luigi.util import requires, inherits
from luigi import LocalTarget

from fieldpathogenomics.luigi.slurm import SlurmExecutableTask
from fieldpathogenomics.utils import CheckTargetNonEmpty
import fieldpathogenomics.utils as utils

import Assemble

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
    lib_list = luigi.ListParameter()
    library = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 1000
        self.n_cpu = len(self.lib_list)
        self.partition = "tgac-medium"

        self.cwd = os.path.join(self.scratch_dir, PIPELINE, VERSION)

    def requires(self):
        return [self.clone(UnzipFastqGZ, library=lib.rstrip()) for lib in self.lib_list]

    def output(self):
        return {lib: [LocalTarget(os.path.join(self.cwd, "nextclip_done", lib + '_nc_ABC_R1.fastq')),
                      LocalTarget(os.path.join(self.cwd, "nextclip_done", lib + '_nc_ABC_R2.fastq'))] for lib in self.lib_list}

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


@inherits(FetchFastqGZ)
@inherits(Assemble.Trimmomatic)
@inherits(LMP_process)
class KatCompPE(CheckTargetNonEmpty, SlurmExecutableTask):

    pe_lib = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 8000
        self.n_cpu = 4
        self.partition = "tgac-medium"

    def requires(self):
        return {'pe': self.clone(Assemble.Trimmomatic, library=self.pe_lib),
                'lmp': self.clone(LMP_process)}

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, self.library, self.library + "-main.mx"))

    def work_script(self):
        return '''#!/bin/bash
                    source kat-2.3.2
                    set -euo pipefail

                    kat comp -n -o {output_prefix} -t {n_cpu} <(cat {lmp_R1} {lmp_R2}) <(zcat {pe_R1} {pe_R2})
                    kat plot spectra-mx -o {output_prefix}.spectra_mx.png -x 100 --intersection {output_prefix}-main.mx


        '''.format(output_prefix=os.path.join(os.path.split(self.output().path)[0], self.library),
                   n_cpu=self.n_cpu,
                   lmp_R1=self.input()['lmp'][self.library][0].path,
                   lmp_R2=self.input()['lmp'][self.library][1].path,
                   pe_R1=self.input()['pe'][0].path,
                   pe_R2=self.input()['pe'][1].path)


@inherits(FetchFastqGZ)
@inherits(Assemble.W2RapContigger)
@inherits(LMP_process)
class MapContigs(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 750
        self.n_cpu = 8
        self.partition = "tgac-medium"

    def requires(self):
        return {'contigs': self.clone(Assemble.W2RapContigger, lib_list=['LIB17363', 'LIB26234']),
                'lmp': self.clone(LMP_process)}

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
                   R1=self.input()['lmp'][self.library][0].path,
                   R2=self.input()['lmp'][self.library][1].path,
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
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, self.library, self.library + ".bam"))

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


@inherits(RawFastQC)
@inherits(KatCompPE)
@inherits(CollectISMetrics)
class PerLibPipeline(luigi.WrapperTask):
    '''Wrapper task that runs all tasks on a single library'''

    def requires(self):
        return [self.clone(RawFastQC),
                self.clone(KatCompPE),
                self.clone(CollectISMetrics)]

    def output(self):
        return self.input()


@inherits(PerLibPipeline)
class LibraryBatchWrapper(luigi.WrapperTask):
    '''Wrapper task to execute the per library part of the pipeline on all
        libraries in :param list lib_list:'''
    library = None

    def requires(self):
        return [self.clone(PerLibPipeline, library=lib.rstrip()) for lib in self.lib_list]

    def output(self):
        return self.input()


@inherits(Assemble.W2RapContigger)
class SOAPPrep(CheckTargetNonEmpty, SlurmExecutableTask):

    prefix = luigi.Parameter(default='pst')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 28000
        self.n_cpu = 1
        self.partition = "tgac-short"

    def requires(self):
        return self.clone(Assemble.W2RapContigger, lib_list=['LIB17363', 'LIB26234'])

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, PIPELINE, VERSION, 'SOAP', self.prefix + '.contig'))

    def work_script(self):
        return '''#!/bin/bash
                    set -euo pipefail
                    cd {cwd}
                    export soap='/usr/users/ga004/buntingd/w2rap/deps/soap_scaffolder'

                    $soap/s_prepare -g {prefix} -K 71 -c {contigs}

        '''.format(contigs=os.path.join(self.input().path, 'a.lines.fasta'),
                   cwd=os.path.split(self.output().path)[0],
                   prefix=self.prefix)


@inherits(LMP_process)
@inherits(Assemble.Trimmomatic)
class SOAPConfig(CheckTargetNonEmpty, luigi.Task):

    def requires(self):
        return {'lmp': self.clone(LMP_process),
                'pe': {lib: self.clone(Assemble.Trimmomatic, library=lib) for lib in ['LIB17363', 'LIB26234']}}

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
        return LocalTarget(os.path.join(self.scratch_dir, PIPELINE, VERSION, 'SOAP', self.prefix + '.readOnContig.gz'))

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
                   prefix=self.prefix)


@requires(SOAPMap)
class SOAPScaff(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 500
        self.n_cpu = 12
        self.partition = "tgac-medium"

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, PIPELINE, VERSION, 'SOAP', self.prefix + '.scaf'))

    def work_script(self):
        return '''#!/bin/bash
                    set -euo pipefail

                    cd {cwd}
                    export soap='/usr/users/ga004/buntingd/w2rap/deps/soap_scaffolder'

                    $soap/s_scaff -p {n_cpu} -g {prefix}

        '''.format(cwd=os.path.split(self.output().path)[0],
                   n_cpu=self.n_cpu,
                   prefix=self.prefix)


@inherits(SOAPScaff)
@inherits(LibraryBatchWrapper)
class Wrapper(luigi.WrapperTask):
    '''Wrapper task that runs all tasks on a single library'''
    library = None

    def requires(self):
        return [self.clone(LibraryBatchWrapper),
                self.clone(SOAPScaff)]

    def output(self):
        return self.input()


if __name__ == '__main__':
    os.environ['TMPDIR'] = "/tgac/scratch/buntingd"
    logger, alloc_log = utils.logging_init(log_dir=os.path.join(os.getcwd(), 'logs'),
                                           pipeline_name=PIPELINE)

    with open(sys.argv[1], 'r') as libs_file:
        lib_list = [line.rstrip() for line in libs_file if line[0] != '#']

    luigi.run(['Wrapper',
               '--lib-list', json.dumps(lib_list),
               '--pe-lib', 'LIB17363'] + sys.argv[2:])
