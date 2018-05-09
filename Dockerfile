FROM karrlab/build:0.0.23

RUN rm ~/.gitconfig

RUN git clone https://github.com/KarrLab/kinetic_datanator.git home/
RUN pip3 install -U --process-dependency-links -e home/[all]
