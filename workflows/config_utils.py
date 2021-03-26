import os
from abc import ABC, abstractmethod

class Data(ABC):
    """
    Abstract class that describes a simple interface for how local/remote
    data will be obtained in snakemake rules.
    """
    def __init__(self, input, output):
        self.input = input
        self.output = output

    @abstractmethod
    def get_cmd(self) -> str:
        """
        returns a shell directive action (eg. downloading, or creating symlink)
        """
        pass


class LocalData(Data):
    def __init__(self, input, output):
        super().__init__(input, output)

    def get_cmd(self) -> str:
        """
        returns symlink command to be run in snakemake shell directive
        """
        return f"ln -sr {self.input} {self.output}"


class S3Data(Data):
    def __init__(self, input, output):
        super().__init__(input, output)

    def get_cmd(self) -> str:
        """
        returns aws s3 copy command to be run in snakemake shell directive
        """
        return f"aws s3 cp {self.input} {self.output}"


def data_factory(config, data_type):
    """
    Create a Data object whose behavior is determined by the config.
    INPUTS:
        * config - snakemake config dictionary
        * data_type - type of data file (eg 'fasta' or 'vcf')
    """
    data_source = config[data_type]['data_source']
    input_file = config[data_type]['file']
    outdir = config['outdir']
    output_file = f'{outdir}/{os.path.basename(input_file)}'

    if data_source == 's3':
        return S3Data(input_file, output_file)
    elif data_source == 'local':
        return LocalData(input_file, output_file)
    else:
        raise ValueError(
            f'Unknown data_source: "{data_source}".'
            ' data_source must be either "s3" or "local".')

        
class Conf:
    """
    Use the snakemake config dict to set variables and configure the Data
    objects which will control the behavior of data aquisition rules.
    """
    def __init__(self, config):
        self.samples = config['samples'].keys() # list of sample names
        self.alignments = config['samples'] # dict of bam/cram indexed by sample
        self.outdir = config['outdir']
        self.delimiter = config['image_filename_delimiter']
        self.fasta = data_factory(config, 'fasta')
        self.fai = data_factory(config, 'fai')
        self.vcf = data_factory(config, 'vcf')

