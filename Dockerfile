FROM ubuntu
FROM perl:latest

# ...
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update
RUN apt-get -y install gcc g++ mono-mcs cmake git subversion gnuplot

ADD . /home

WORKDIR /home

RUN git clone git://git.code.sf.net/p/scalpel/code scalpel
WORKDIR /home/scalpel
RUN make

WORKDIR /usr/local/bin
RUN mkdir ./bamtools-2.3.0/
WORKDIR /usr/local/bin/bamtools-2.3.0/
RUN mkdir ./bin/
WORKDIR /home/scalpel
RUN tar -zxvf protocol_bundle-0.5.3.tar.gz

RUN ln -s /home/scalpel/bamtools-2.3.0/bin/bamtools /usr/local/bin/bamtools-2.3.0/bin
RUN ln -s /home/scalpel/scalpel-discovery /usr/local/bin
RUN ln -s /home/scalpel/scalpel-export /usr/local/bin
RUN ln -s /home/scalpel/FindSomatic.pl /usr/local/bin
RUN ln -s /home/scalpel/FindDenovos.pl /usr/local/bin
RUN ln -s /home/scalpel/FindVariants.pl /usr/local/bin
RUN ln -s /home/scalpel/HashesIO.pm /usr/local/bin
RUN ln -s /home/scalpel/Usage.pm /usr/local/bin
RUN ln -s /home/scalpel/Utils.pm /usr/local/bin
RUN ln -s /home/scalpel/Microassembler /usr/local/bin
RUN ln -s /home/scalpel/protocol_bundle-0.5.3 /usr/local/bin

#----------------------------------------------------------------------------------

ENV HELP='\nusage:\n\tscalpel-discovery <COMMAND> [OPTIONS]\n\tscalpel-export <COMMAND> [OPTIONS]'

CMD echo $HELP