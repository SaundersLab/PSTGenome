import os
import sys
import json
import tempfile

import luigi

from luigi.util import requires, inherits
from luigi import LocalTarget
from luigi.file import TemporaryFile


from fieldpathogenomics.luigi.slurm import SlurmExecutableTask
from fieldpathogenomics.utils import CheckTargetNonEmpty
import fieldpathogenomics.utils as utils

import Assemble

PIPELINE = os.path.basename(__file__).split('.')[0]
VERSION = '0.1'
luigi.auto_namespace(scope=PIPELINE)


# This information has to be filled in based on the output of
# CollectISMetrics, obviously hardcoding it here is not good!!!

INSERT_SIZES = {'LIB17363': 400, 'LIB26234': 800,
                'LIB19826': 0, 'LIB19827': 0,
                'LIB19828': 0, 'LIB19829': 0,
                'LIB19830': 0, 'LIB19831': 0,
                'LIB19832': 6000, 'LIB19833': 5000,
                'LIB19834': 4000, 'LIB19835': 3500,
                'LIB19836': 2500, 'LIB19837': 2100}

# ------------------ Homozygous filter -------------------------- #


@requires(Assemble.W2RapContigger)
class Redundans(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 4000
        self.n_cpu = 4
        self.partition = "tgac-medium"

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
                              --overlap 0.95  \
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
class BWAIndexContigs(Assemble.BWAIndex):
    pass


# ------------------ Scaffolding -------------------------- #


@requires(HomozygousContigs)
class SOAPPrep(CheckTargetNonEmpty, SlurmExecutableTask):

    prefix = luigi.Parameter(default='PST')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 28000
        self.n_cpu = 1
        self.partition = "tgac-short"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, 'SOAP', 'K' + str(self.K), self.prefix + str(self.K) + '.contig'))

    def work_script(self):
        return '''#!/bin/bash
                    set -euo pipefail
                    cd {cwd}
                    export soap='/usr/users/ga004/buntingd/w2rap/deps/soap_scaffolder'

                    $soap/s_prepare -g {prefix} -K 71 -c {contigs}

        '''.format(contigs=self.input().path,
                   cwd=os.path.dirname(self.output().path),
                   prefix=self.prefix + str(self.K))


@inherits(Assemble.LMP_process)
@inherits(Assemble.Trimmomatic)
class SOAPConfig(CheckTargetNonEmpty, luigi.Task):

    def requires(self):
        return {'lmp': self.clone(Assemble.LMP_process),
                'pe': {lib: self.clone(Assemble.Trimmomatic, library=lib) for lib in self.pe_libs}}

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
        self.mem = 2000
        self.n_cpu = 12
        self.partition = "tgac-medium"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, 'SOAP', 'K' + str(self.K), self.prefix + str(self.K) + '.readOnContig.gz'))

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
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, 'SOAP', 'K' + str(self.K), self.prefix + str(self.K) + '.scafSeq'))

    def work_script(self):
        return '''#!/bin/bash
                    set -euo pipefail

                    cd {cwd}
                    export soap='/usr/users/ga004/buntingd/w2rap/deps/soap_scaffolder'

                    $soap/s_scaff -p {n_cpu} -g {prefix}

        '''.format(cwd=os.path.split(self.output().path)[0],
                   n_cpu=self.n_cpu,
                   prefix=self.prefix + str(self.K))


@requires(SOAPScaff)
class RedundansScaffoldStats(Assemble.AbyssFac):
    pass


@requires(SOAPScaff)
class SOAPNremap(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 6000
        self.n_cpu = 1
        self.partition = "tgac-medium"
        self.rm_tmp = False

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, 'SOAP', 'K' + str(self.K), self.prefix + str(self.K) + '.scaf.fasta'))

    def work_script(self):
        soap_base = self.input().path[:-8]

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

# ------------------ Gap Filling -------------------------- #


@inherits(SOAPNremap)
@inherits(Assemble.AbyssBloomBuild)
class AbyssSealerReduced(CheckTargetNonEmpty, SlurmExecutableTask):

    sealer_klist = luigi.ListParameter()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 2500
        self.n_cpu = 5
        self.partition = "tgac-medium"

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, "sealer", "SOAP", 'K' + str(self.K), 'K' + str(self.K) + "_scaffold.fa"))

    def requires(self):
        return {'bloomfilters': [self.clone(Assemble.AbyssBloomBuild, bloom_k=k) for k in self.sealer_klist],
                'scaffolds': self.clone(SOAPNremap)}

    def work_script(self):
        return '''#!/bin/bash
                    source abyss-2.0.2;
                    mkdir -p {output}/temp
                    set -euo pipefail

                    abyss-sealer {k_args} -P25 --flank-length=150 -j {n_cpu} -o {output}/temp/{prefix} -S {scaffolds} {bloomfilters}

                    mv {output}/temp/{prefix}* {output}/
        '''.format(k_args=' '.join(['-k' + str(k) for k in self.sealer_klist]),
                   bloomfilters=' '.join(['-i ' + x.path for x in self.input()['bloomfilters']]),
                   scaffolds=self.input()['scaffolds'].path,
                   n_cpu=self.n_cpu,
                   output=os.path.dirname(self.output().path),
                   prefix='K' + str(self.K))


@requires(AbyssSealerReduced)
class BWAIndexScaffolds(Assemble.BWAIndex):
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
        self.partition = "tgac-medium"

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


@inherits(BWAIndexScaffolds)
@inherits(BAMtoFASTQ)
class MapLinkedReadsScaffolds(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 500
        self.n_cpu = 24
        self.partition = "tgac-medium"

    def requires(self):
        return {'contigs': self.clone(BWAIndexScaffolds),
                'reads': self.clone(BAMtoFASTQ)}

    def output(self):
        return LocalTarget(os.path.join(self.scratch_dir, PIPELINE, VERSION, "linked_reads", 'scaffolds', "K" + str(self.K), "barcoded_aligned.bam"))

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


@inherits(MapLinkedReadsScaffolds)
@inherits(SOAPNremap)
class ArcsLinksScaffolds(CheckTargetNonEmpty, SlurmExecutableTask):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set the SLURM request params for this task
        self.mem = 8000
        self.n_cpu = 1
        self.partition = "tgac-medium"

    def requires(self):
        return {'contigs': self.clone(SOAPNremap),
                'reads': self.clone(MapLinkedReadsScaffolds)}

    def output(self):
        return LocalTarget(os.path.join(self.base_dir, PIPELINE, VERSION, 'linked_reads', 'scaffolds', "K" + str(self.K), "K" + str(self.K) + ".scaffolds.fa"))

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
class RedundansARCSStats(Assemble.AbyssFac):
    pass


@inherits(RedundansScaffoldStats)
@inherits(RedundansARCSStats)
@inherits(AbyssSealerReduced)
class Wrapper(luigi.WrapperTask):
    library = None
    K = None
    bloom_k = None
    bloom_size = None
    K_list = luigi.ListParameter(default=[280])

    def requires(self):
        for k in self.K_list:
            yield self.clone(RedundansScaffoldStats, K=k)
            yield self.clone(RedundansARCSStats, K=k)
            yield self.clone(AbyssSealerReduced, K=k, sealer_klist=self.sealer_klist)

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
               '--sealer-klist', json.dumps([30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 128]),
               '--pe-libs', json.dumps(pe_libs),
               '--lmp-libs', json.dumps(lmp_libs)] + sys.argv[3:])
