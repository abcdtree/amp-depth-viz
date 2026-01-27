# amp-depth-viz
## visualise amplicon genome coverage by bokeh using mosdepth output

## install
```
    conda create -n amp-depth-viz -c conda-forge -c bioconda -c nanoporetech fastcat python=3.14
    conda activate amp-depth-viz
    pip install amp-depth-viz
```

## Instruction
```
    #running fastcat and provide a samplesheet
    #fastq_pass is the fastq_pass folder from the sars ONT run
    #samplesheet requests at least two columns: sampleID,barcode
    #It helps to filter the barcode and relate to right sampleID
    amp-depth-viz sample/sample_input.bed --fastq_pass Data/sars/fastq_pass --samplesheet samplesheet.csv --output test.html

    #if you already ran fastcat and have the per reads summary file, you can do 
    amp-depth-viz sample/sample_input.bed --fastcat_perreads Data/per-reads-summary.tsv --samplesheet samplesheet.csv --output test.html
```

## Using import inner functions
```
    #import functions inside amp-depth-viz package
    from amp_depth_viz.main import *

    #load template from the package files -- add in 0.1.5
    template = load_template()

    #create plots
    fig_len, fig_qual, fig_cnt = plot_summary_plotly(file_path_to_fastcat_per_read_file)
    coverage_plot_plotly = plot_coverage_plotly(
            coverage_bed_file_from_mosdepth,
            threshold=20,
            xlim=30000,
            ylim=800,
            ncols=3,
            colors=Colors(),
        )
      
    #combine the plots and create html
    plotly_plot_combine([fig_len, fig_qual, fig_cnt, coverage_plot_plotly], template=template, output_path=output_path_to_save_html)
```

## Usage
```
    usage: amp-depth-viz [-h] [--fastcat_perreads FASTCAT_PERREADS | --fastq_pass FASTQ_PASS] [--template TEMPLATE]
                     [--samplesheet SAMPLESHEET] [--output OUTPUT] [--xlim XLIM] [--ylim YLIM] [--threshold THRESHOLD]
                     [--threads THREADS] [--ncols NCOLS]
                     coveragebed

create wf-artic like amplicon coverage plots using bokeh

positional arguments:
  coveragebed           a full sample bed files contains all genome depth information from mosdepth

options:
  -h, --help            show this help message and exit
  --fastcat_perreads FASTCAT_PERREADS
                        per-read-stats.tsv output from fastcat
  --fastq_pass FASTQ_PASS
                        the fastq_folder from a sars cov run
  --template TEMPLATE   jinja2 template to render -- default in template/template.html
  --samplesheet SAMPLESHEET
                        An optional samplesheet to rename the barcode, at least two columns required: sampleID, barcode
  --output OUTPUT       html output file
  --xlim XLIM           max of the genome position
  --ylim YLIM           the expected max depth
  --threshold THRESHOLD
                        depth threshold for passing QC , default as 20
  --threads THREADS     max number of cpus to use for fastcat analysis
  --ncols NCOLS         number of columns to grid the plots in the html pages, default as 3
  --plot_package, -p {bokeh,plotly}
                        Choose a python plot package to make the plots,from the allowed options: bokeh, plotly
```