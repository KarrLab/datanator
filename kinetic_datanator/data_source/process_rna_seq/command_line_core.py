from kinetic_datanator.data_source.process_rna_seq import core
import cement

class BaseController(cement.Controller):

    class Meta:
        label = 'base'
        description = 'Utilities for aggregating data for biochemical models'


class DownloadCDNA(cement.Controller):

    class Meta:
        label = 'download-cdna'
        description = "Download CDNA file into specified directory"
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['strain_name'], dict(type=str, help="strain of the reference genome")),
            (['url'], dict(type=str, help="url to download the reference genome")),
            (['temp_directory'], dict(type=str, help="location to save the reference genome")),
        ]

    @cement.ex(hide=True)
    def _default(self):
        core.download_cdna(self.app.pargs.strain_name, self.app.pargs.url, self.app.pargs.temp_directory)


class DownloadFASTQ(cement.Controller):

    class Meta:
        label = 'download-fastq'
        description = "Download FASTQ files into specified directory"
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['experiment_name'], dict(type=str, help="Name of the experiment")),
            (['sample_name'], dict(type=str, help="Name of the sample")),
            (['temp_directory'], dict(type=str, help="location to save the fastq files")),
            (['fastq_urls'], dict(type=str, help="space seperated list of fastq urls")),
        ]
    @cement.ex(hide=True)
    def _default(self):
        pargs = self.app.pargs
        core.download_fastq(pargs.experiment_name, pargs.sample_name, pargs.temp_directory, pargs.fastq_urls)



class ProcessCDNA(cement.Controller):

    class Meta:
        label = 'process-cdna'
        description = "Process a CDNA file into a kallisto index file and save it to a specified directory"
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['strain_name'], dict(type=str, help="strain of the reference genome")),
            (['output_directory'], dict(type=str, help="location to save index file")),
            (['temp_directory'], dict(type=str, help="location of the reference genomes")),
        ]

    @cement.ex(hide=True)
    def _default(self):
        pargs = self.app.pargs
        core.process_cdna(pargs.strain_name, pargs.output_directory, pargs.temp_directory)


class ProcessFASTQ(cement.Controller):

    class Meta:
        label = 'process-fastq'
        description = "Process FASTQ files"
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['experiment_name'], dict(type=str, help="Name of the experiment")),
            (['sample_name'], dict(type=str, help="Name of the sample")),
            (['strain_name'], dict(type=str, help="strain name of the reference genome")),
            (['num_fastq_files'], dict(type=int, help="Number of FASTQ files for sample")),
            (['read_type'], dict(type=str, help="Read type of RNA-seq experiment")),
            (['output_directory'], dict(type=str, help="location to save the processed files")),
            (['temp_directory'], dict(type=str, help="location of the fastq files")),
        ]
    @cement.ex(hide=True)
    def _default(self):
        pargs = self.app.pargs
        core.process_fastq(pargs.experiment_name, pargs.sample_name, pargs.strain_name, pargs.num_fastq_files, pargs.read_type, pargs.output_directory, pargs.temp_directory)

class DeleteCDNAFiles(cement.Controller):

    class Meta:
        label = 'delete-cdna'
        description = "Delete CDNA files"
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['strain_name'], dict(type=str, help="Name of the strain")),
            (['temp_directory'], dict(type=str, help="location where strain is saved")),

        ]
    @cement.ex(hide=True)
    def _default(self):
        pargs = self.app.pargs
        core.delete_cdna_files(pargs.strain_name, pargs.temp_directory)

class DeleteFASTQFiles(cement.Controller):

    class Meta:
        label = 'delete-fastq'
        description = "Delete fastq files"
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['experiment_name'], dict(type=str, help="Name of the experiment")),
            (['sample_name'], dict(type=str, help="Name of the sample")),
            (['temp_directory'], dict(type=str, help="location where strain is saved")),

        ]
    @cement.ex(hide=True)
    def _default(self):
        pargs = self.app.pargs
        core.delete_fastq_files(pargs.experiment_name, pargs.sample_name, pargs.temp_directory)


class App(cement.App):

    class Meta:
        label = "kinetic_datanator"
        base_controller = "base"
        handlers = [
            BaseController,
            DownloadCDNA,
            DownloadFASTQ,
            ProcessCDNA,
            ProcessFASTQ,
            DeleteCDNAFiles,
            DeleteFASTQFiles,

        ]

def main():
    with App() as app:
        app.run()

if __name__ == "__main__":
    main()

