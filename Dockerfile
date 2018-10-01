FROM karrlab/wc_env_dependencies:latest

<<<<<<< HEAD
ADD . /home
WORKDIR /home

RUN pip3 install -U --process-dependency-links -e /home/[all]


CMD python3 manage.py runserver --host 0.0.0.0 --port ${PORT}
=======
RUN git clone https://github.com/KarrLab/datanator.git home/
RUN pip3 install -U --process-dependency-links -e home/[all]
>>>>>>> 650d74fcc6d11695e9a12ba7db13ea71ad681d51
