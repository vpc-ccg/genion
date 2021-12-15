FROM ubuntu:16.04

RUN apt-get update -y --fix-missing
RUN apt-get install git make g++ zlib1g-dev -y

RUN mkdir /genion
RUN mkdir /input
RUN mkdir /output
WORKDIR /genion

RUN git clone https://github.com/vpc-ccg/genion.git /genion --recursive -b docker
RUN make

ENV PATH="/genion:${PATH}"
VOLUME /input
VOLUME /output
ENTRYPOINT ["/genion/genion"]
