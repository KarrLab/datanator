FROM lzy7071/karrlabdatanator_dependencies:latest

COPY . /home
WORKDIR /home

RUN pip3 install -e .

ENTRYPOINT ["python3"]
CMD ["/home/datanator/rest/__init__.py"]