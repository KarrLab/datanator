FROM lzy7071/karrlabdatanator_dependencies:latest

RUN mkdir -p /tmp/datanator

RUN apt-get update -y \
	&& apt-get install -y --no-install-recommends \
		curl \
		wget \
	&& rm -rf /var/lib/apt/lists/*

WORKDIR /root
ENTRYPOINT ["python3"]
CMD ["/home/datanator/rest/__init__.py"]
