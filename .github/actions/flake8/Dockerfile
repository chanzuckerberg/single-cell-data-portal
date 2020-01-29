FROM python:3

LABEL "com.github.actions.name"="Python linting and code q"
LABEL "com.github.actions.description"="Python linting using black and flake8"
LABEL "com.github.actions.icon"="bell"
LABEL "com.github.actions.color"="purple"

COPY requirements.txt /requirements.txt
COPY .flake8 /.flake8
RUN pip3 install -r /requirements.txt

COPY entrypoint.sh /entrypoint.sh
ENTRYPOINT ["/entrypoint.sh"]