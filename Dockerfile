FROM lzy7071/karrlabdatanator_dependencies:latest

RUN mkdir -p /tmp/datanator

COPY . /tmp/datanator

RUN cd /tmp/datanator \
	&& pip3 install -e . \
	&& rm -rf /tmp/datanator

WORKDIR /root
ENTRYPOINT ["python3"]
CMD ["/home/datanator/rest/__init__.py"]
