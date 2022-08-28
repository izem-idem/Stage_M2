FROM rocker/tidyverse:4.1.3

#libcurl4-gnutils-dev si j'en ai besoin jle resort
RUN apt-get update --allow-releaseinfo-change && apt-get install -y python3 python3-pip  python3-tk libpq-dev build-essential libxml2-dev libssl-dev


#je ne sais pas prk je dois mettre mon script d'installation dans /opt/software/setup/R
RUN mkdir -p /opt/software/setup/R
ADD install.R /opt/software/setup/R/
RUN Rscript /opt/software/setup/R/install.R

RUN mkdir /app

WORKDIR /app

COPY requirements.txt .
RUN pip3 install -r requirements.txt
#copy Explore dir
ADD . /app

CMD python3 Main.py


