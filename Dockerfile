FROM karrlab/wc_env_dependencies:latest

ADD . /home
WORKDIR /home

RUN pip3 install -U --process-dependency-links -e /home/[all]


CMD python3 manage.py runserver --host 0.0.0.0 --port ${PORT}
