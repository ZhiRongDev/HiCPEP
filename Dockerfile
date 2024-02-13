FROM jupyter/base-notebook:python-3.11.5

USER root

COPY . /tmp
WORKDIR /tmp

RUN \
	DEBIAN_FRONTEND=noninteractive apt-get update && \
	DEBIAN_FRONTEND=noninteractive apt-get install -y build-essential

RUN pip install -r requirements.txt

WORKDIR "${HOME}/work"

USER 1001
