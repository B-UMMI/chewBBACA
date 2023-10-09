FROM ubuntu:22.04
MAINTAINER UMMI, imm-bioinfo@medicina.ulisboa.pt

# Create directory to store files
WORKDIR /WORKY/

# Install curl and wget to download files, git to clone repos and Python 3 to install and use chewBBACA
RUN apt update && apt install -y wget curl git python3 python3-setuptools

# Install BLAST
WORKDIR /WORKY/BLAST
# Download and extract
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz && \
	tar -xzvf ncbi-blast-2.9.0+-x64-linux.tar.gz && \
	rm ncbi-blast-2.9.0+-x64-linux.tar.gz
# Add BLAST executables to PATH
ENV PATH="/WORKY/BLAST/ncbi-blast-2.9.0+/bin:${PATH}"

# Install MAFFT
WORKDIR /WORKY/MAFFT
RUN wget https://mafft.cbrc.jp/alignment/software/mafft_7.520-1_amd64.deb
# Install
RUN dpkg -i mafft_7.520-1_amd64.deb
ENV PATH="/WORKY/MAFFT:${PATH}"

# Install FastTree
WORKDIR /WORKY/FastTree
RUN wget http://www.microbesonline.org/fasttree/FastTree
# Change permissions for executable
RUN chmod 755 *
ENV PATH="/WORKY/FastTree:${PATH}"

# INSTALL chewBBACA and dependencies listed in requirements.txt
WORKDIR /WORKY/chewBBACA
# Download latest chewBBACA release
RUN curl -ks https://api.github.com/repos/B-UMMI/chewBBACA/releases/latest | \
	grep "tarball_url" | \
	cut -d "\"" -f4 | \
	xargs -I {} wget --no-check-certificate -O chewBBACA.tar.gz {}
RUN tar -xf chewBBACA.tar.gz
RUN mv B-UMMI-chewBBACA*/* .
RUN rm -rf B-UMMI-chewBBACA*
# Install chewBBACA and dependencies
RUN python3 setup.py install
ENV PATH="/WORKY/chewBBACA/chewBBACA:${PATH}"

WORKDIR /WORKY/
