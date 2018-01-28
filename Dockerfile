FROM karrlab/build:0.0.6

RUN pip3 install -U -e git+git://github.com/KarrLab/wc_utils#egg=wc_utils
RUN pip3 install unicodecsv
RUN pip3 install pandas

ADD . /home/
