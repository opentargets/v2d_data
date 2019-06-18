FROM continuumio/miniconda:4.6.14

COPY ./environment.yaml /v2d/
WORKDIR /v2d
RUN conda env create -n v2d_data --file environment.yaml
RUN echo "source activate v2d_data" > ~/.bashrc
ENV PATH /opt/conda/envs/v2d_data/bin:$PATH

RUN curl https://sdk.cloud.google.com | bash && echo "source /root/google-cloud-sdk/completion.bash.inc; source /root/google-cloud-sdk/path.bash.inc" >> ~/.bashrc

COPY ./ /v2d
