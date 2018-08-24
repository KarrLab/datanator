FROM karrlab/wc_env_dependencies:latest

RUN git clone https://github.com/KarrLab/kinetic_datanator.git home/
RUN pip3 install -U --process-dependency-links -e home/[all]
