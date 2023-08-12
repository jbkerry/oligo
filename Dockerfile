FROM python:3.9 as builder

WORKDIR /usr/local

COPY docker/RepeatMasker-4.1.5.tar.gz .

RUN tar -xvzf RepeatMasker-4.1.5.tar.gz --exclude=RepeatMasker/Libraries/Dfam.h5

COPY docker/Dfam.h5 /usr/local/RepeatMasker/Libraries/Dfam.h5

RUN apt-get update && apt-get install -y rsync

RUN wget https://github.com/Benson-Genomics-Lab/TRF/archive/refs/tags/v4.09.1.tar.gz \
    && tar -xvzf v4.09.1.tar.gz \
    && cd TRF-4.09.1 \
    && mkdir build \
    && cd build \
    && ../configure \
    && make \
    && rm ../../v4.09.1.tar.gz

RUN mkdir blat \
    && cd blat \
    && rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/blat/blat ./

FROM python:3.9-slim

WORKDIR /usr/local

COPY --from=builder /usr/local/RepeatMasker ./RepeatMasker
COPY --from=builder /usr/local/blat/blat ./blat
COPY --from=builder /usr/local/TRF-4.09.1 ./TRF-4.09.1

COPY docker/humans.hmm .
COPY docker/mouse.hmm .
COPY config.txt .

RUN apt-get update && apt-get install -y hmmer build-essential libcurl4

ENV PERL5LIB="$PERL5LIB:/usr/local/RepeatMasker"

ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
RUN pip install --upgrade pip \
    && pip install oligo-capture

RUN cd RepeatMasker \
    && perl configure -trf_prgm /usr/local/TRF-4.09.1/build/src/trf -hmmer_dir /usr/bin -default_search_engine hmmer

RUN mkdir /results
WORKDIR /results

ENTRYPOINT ["python", "-m", "oligo", "-cfg", "/usr/local/config.txt"]