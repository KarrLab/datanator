FROM lzy7071/karrlabdatanator_dependencies:latest

ADD . /home
WORKDIR /home

RUN pip3 install -U -e /home/[all]

ENTRYPOINT ["python3"]
CMD ["/home/datanator/rest/__init__.py"]