FROM python:3.9

# change root password
RUN echo 'root:root' | chpasswd

# install OGSync
ADD OGSync /opt/OGSync/
ADD oglib/*.py /opt/OGSync/oglib/

# add PATH
ENV PATH $PATH:/opt/OGSync
ENV PYTHONPATH="${PYTHONPATH}:/usr/local/lib/python3.9/site-packages"

# install Python modules
RUN pip3 install --upgrade pip
RUN pip3 install colorama pymongo requests progressbar biopython prettytable

# init command
CMD ["/bin/bash"]