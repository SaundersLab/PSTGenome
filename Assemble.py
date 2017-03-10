import os
import sys
import json

import luigi
from luigi.util import requires, inherits
from luigi import LocalTarget

from fieldpathogenomics.luigi.slurm import SlurmExecutableTask
from fieldpathogenomics.utils import CheckTargetNonEmpty
import fieldpathogenomics.utils as utils


PIPELINE = os.path.basename(__file__).split('.')[0]
VERSION = '0.0'


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
                           adapters='/tgac/software/testing/trimmomatic/0.30/x86_64/bin/adapters/TruSeq.cat.fa',
                           R1_out=self.output()[0].path,
                           R2_out=self.output()[1].path)


@requires(FetchFastqGZ)
class FastxQC(SlurmExecutableTask):
    '''Runs Fastx toolkit to plot the nucleotide and base call quality score distributions '''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 2000
        self.n_cpu = 1
        self.partition = "tgac-medium"

    def output(self):
        working_dir = os.path.join(self.base_dir, PIPELINE, VERSION, self.library)
        return {'stats_R1': LocalTarget(os.path.join(working_dir, 'QC', self.library + "_R1_stats.txt")),
                'stats_R2': LocalTarget(os.path.join(working_dir, 'QC', self.library + "_R2_stats.txt")),
                'boxplot_R1': LocalTarget(os.path.join(working_dir, 'QC', self.library + "_R1_quality.png")),
                'boxplot_R2': LocalTarget(os.path.join(working_dir, 'QC', self.library + "_R2_quality.png")),
                'nt_dist_R1': LocalTarget(os.path.join(working_dir, 'QC', self.library + "_R1_nt_distr.png")),
                'nt_dist_R2': LocalTarget(os.path.join(working_dir, 'QC', self.library + "_R2_nt_distr.png")),
                }

    def work_script(self):
        return '''#!/bin/bash
        source fastx_toolkit-0.0.13.2
        set -euo pipefail

        gzip -cd {R1_in} | fastx_quality_stats -o {stats_R1} -Q33
        gzip -cd {R2_in} | fastx_quality_stats -o {stats_R2} -Q33

        fastq_quality_boxplot_graph.sh -i {stats_R1} -o {boxplot_R1}
        fastq_quality_boxplot_graph.sh -i {stats_R2} -o {boxplot_R2}

        fastx_nucleotide_distribution_graph.sh -i {stats_R1} -o {nt_dist_R1}
        fastx_nucleotide_distribution_graph.sh -i {stats_R2} -o {nt_dist_R2}

        '''.format(R1_in=self.input()[0].path,
                   R2_in=self.input()[1].path,
                   stats_R1=self.output()['stats_R1'].path,
                   stats_R2=self.output()['stats_R2'].path,
                   boxplot_R1=self.output()['boxplot_R1'].path,
                   boxplot_R2=self.output()['boxplot_R2'].path,
                   nt_dist_R1=self.output()['nt_dist_R1'].path,
                   nt_dist_R2=self.output()['nt_dist_R2'].path)


@requires(Trimmomatic)
class FastQC(CheckTargetNonEmpty, SlurmExecutableTask):

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            # Set the SLURM request params for this task
            self.mem = 2000
            self.n_cpu = 1
            self.partition = "tgac-medium"

        def output(self):
            return [LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, self.library, 'QC', 'R1', 'fastqc_data.txt')),
                    LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, self.library, 'QC', 'R2', 'fastqc_data.txt'))]

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
                    '''.format(output_dir=os.path.join(self.scratch_dir, PIPELINE, VERSION, self.library, 'FastQC'),
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
        self.n_cpu = 8
        self.partition = "tgac-medium"

    def output(self):
        return [LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, self.library + ".hist")),
                LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, self.library + ".hist.png"))]

    def work_script(self):
        return '''#!/bin/bash
                    source kat-2.3.2
                    set -euo pipefail

                    kat hist -o {output_prefix} -t {n_cpu} <(zcat {R1} {R2})

        '''.format(output_prefix=self.library + '.hist',
                   n_cpu=self.n_cpu,
                   R1=self.input()[0].path,
                   R2=self.input()[1].path)


@inherits(KatHist)
@inherits(FastQC)
@inherits(FastxQC)
class PerLibPipeline(luigi.WrapperTask):
    '''Wrapper task that runs all tasks on a single library'''

    def requires(self):
        return {'kathit': self.clone(KatHist),
                'fastqc': self.clone(FastQC),
                'fastxgq': self.clone(FastxQC)}

    def output(self):
        return self.input()


@inherits(PerLibPipeline)
class LibraryBatchWrapper(luigi.WrapperTask):
    '''Wrapper task to execute the per library part of the pipeline on all
        libraries in :param list lib_list:'''
    lib_list = luigi.ListParameter()
    library = None

    def requires(self):
        return [self.clone_parent(library=lib.rstrip()) for lib in self.lib_list]

    def output(self):
        return self.input()

# ----------------------------------------------------------------------- #


if __name__ == '__main__':
    os.environ['TMPDIR'] = "/tgac/scratch/buntingd"
    logger, alloc_log = utils.logging_init(log_dir=os.path.join(os.getcwd(), 'logs'),
                                           pipeline_name=PIPELINE)

    with open(sys.argv[1], 'r') as libs_file:
        lib_list = [line.rstrip() for line in libs_file if line[0] != '#']

    luigi.run(['LibraryBatchWrapper',
               '--lib-list', json.dumps(lib_list)] + sys.argv[2:])
