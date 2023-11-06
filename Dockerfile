FROM jupyter/datascience-notebook:6d91620dbb99

RUN conda install -c bioconda \
    trimmomatic \
    fastqc \
    star

RUN pip install HTSeq

#RUN conda install -c bioconda sra-tools==3.0.7

RUN wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz \
    && tar -vxzf sratoolkit.tar.gz \
    && mv sratoolkit.3.0.7-ubuntu64/ /opt/conda/envs/

RUN conda install -c bioconda bioconductor-deseq2

ENV PATH /opt/conda/envs/trimmomatic/bin:$PATH
ENV PATH /opt/conda/envs/fastqc/bin:$PATH
ENV PATH /opt/conda/envs/sratoolkit.3.0.7-ubuntu64/bin:$PATH
ENV PATH /opt/conda/envs/star/bin:$PATH


