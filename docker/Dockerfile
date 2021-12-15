FROM gcc:11.2 

RUN apt-get update -y --fix-missing
RUN apt-get install git make zlib1g-dev -y


RUN mkdir /input
RUN mkdir /output

RUN git clone https://github.com/vpc-ccg/genion.git /genion --recursive -b docker

WORKDIR /genion
RUN make

ENV PATH="/genion:${PATH}"
VOLUME /input
VOLUME /output
ENTRYPOINT ["/genion/genion"]

