FROM gapsystem/gap-docker-master:francy

RUN cd /home/gap/inst/gap-master/pkg && git clone https://github.com/gap-packages/FrancyMonoids.git

COPY . /home/gap
