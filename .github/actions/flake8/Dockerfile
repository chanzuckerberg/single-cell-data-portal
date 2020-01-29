FROM python:3

COPY requirements.txt /requirements.txt
RUN pip3 install -r /requirements.txt

COPY entrypoint.sh /entrypoint.sh
ENTRYPOINT ["/entrypoint.sh"]