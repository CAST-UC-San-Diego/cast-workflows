FROM ubuntu:20.04

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -qqy \
    python3 \
    python3-pip \
    wget \
    build-essential \
    curl \
    unzip \
    && rm -rf /var/lib/apt/lists/*


# Install Python dependencies
RUN pip3 install pandas

# Download and install gcloud
RUN curl https://dl.google.com/dl/cloudsdk/release/google-cloud-sdk.tar.gz > /tmp/google-cloud-sdk.tar.gz
RUN mkdir -p /usr/local/gcloud \
  && tar -C /usr/local/gcloud -xvf /tmp/google-cloud-sdk.tar.gz \
  && /usr/local/gcloud/google-cloud-sdk/install.sh
ENV PATH $PATH:/usr/local/gcloud/google-cloud-sdk/bin

# Copy python script 
COPY convert_phenotype_plink.py /usr/bin/
RUN chmod +x /usr/bin/convert_phenotype_plink.py

# Download and install PLINK2
RUN wget https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20241124.zip \
    && unzip plink2_linux_x86_64_20241124.zip\
    && mv plink2 /usr/bin/ \
    && chmod +x /usr/bin/plink2 \
    && rm plink2_linux_x86_64_20241124.zip









