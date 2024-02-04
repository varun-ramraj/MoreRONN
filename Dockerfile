# Pull base image
FROM python:3-slim-bullseye

# Set environment variables
ENV PIP_DISABLE_PIP_VERSION_CHECK 1
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

# Set work directory
WORKDIR /moreronn_app

# Copy project
COPY . .

# Compile MoreRONN
# and copy-symlink it
RUN apt-get update && apt-get -y install build-essential libtool autoconf automake && ./bootstrap.sh && ./configure && make && make install && cp /usr/local/bin/moreRONN_49 web/bin/moreRONN && cp /moreronn_app/data/*.dat web/bin/
