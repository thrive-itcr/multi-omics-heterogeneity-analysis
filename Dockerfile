# Copyright (c) General Electric Company, 2019.  All rights reserved.

FROM thriveitcr/rt106-algorithm-sdk

USER root

# Setup Java
RUN apt-get -y update \
    && apt-get install -y default-jre

# Setup Python MySQL connection
RUN pip install mysql-connector

RUN mkdir -p /rt106/data
ADD MOHAtool.jar rt106SpecificAdaptorCode.py rt106SpecificAdaptorDefinitions.json /rt106/

RUN chmod a+x /rt106/entrypoint.sh

WORKDIR /rt106
RUN chown -R rt106:rt106 /rt106
USER rt106:rt106

CMD ["/rt106/entrypoint.sh"]
