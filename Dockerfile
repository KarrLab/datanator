FROM karrlab/wc_env_dependencies:latest

ADD . /home
WORKDIR /home

RUN pip3 install -U --process-dependency-links -e /home/[all]


CMD gunicorn -w 4 -b 0.0.0.0:${PORT} --timeout 120 manage:app
