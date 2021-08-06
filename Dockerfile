ARG CORE_VERSION=dev

FROM askcos/askcos-core:$CORE_VERSION as core

USER root

COPY requirements.txt requirements.txt
RUN pip install --no-cache-dir -r requirements.txt && rm requirements.txt

COPY --chown=askcos:askcos . /usr/local/askcos-site

USER askcos

ENV PYTHONPATH=/usr/local/askcos-site:${PYTHONPATH}

RUN python /usr/local/askcos-site/manage.py collectstatic --noinput

LABEL site.version={VERSION} \
      site.git.hash={GIT_HASH} \
      site.git.date={GIT_DATE} \
      site.git.describe={GIT_DESCRIBE}
