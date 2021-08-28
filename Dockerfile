FROM continuumio/miniconda3:4.9.2

RUN apt-get install -y \
    libtbb2 \
 && conda install --channel bioconda \
    cd-hit=4.8.1 \
    hmmer=3.3.2 \
    bowtie2=2.4.1 \
    samtools=1.11 \
    bedtools=2.30.0 \
 && conda clean --all --yes \
 && pip install --no-cache-dir \
    numpy==1.19.2 \
    pandas==1.2.4 \
    matplotlib==3.3.2 \
    seaborn==0.11.0 \
    scipy==1.5.2 \
    statsmodels==0.12.0 \
    argparse==1.1

RUN wget http://cab.spbu.ru/files/release3.15.3/SPAdes-3.15.3-Linux.tar.gz \
 && tar -xzf SPAdes-3.15.3-Linux.tar.gz

ENV PATH="/SPAdes-3.15.3-Linux/bin:${PATH}"

COPY ./MetaGPA/* /MetaGPA/MetaGPA/
COPY ./__main__.py /MetaGPA/
WORKDIR /
ENTRYPOINT ["python", "MetaGPA"]
