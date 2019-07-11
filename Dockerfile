FROM lzy7071/karrlabdatanator_dependencies:latest

RUN mkdir -p /tmp/datanator

COPY . /tmp/datanator

RUN cd /tmp/datanator \
	&& pip3 install -e .

WORKDIR /root
CMD bash