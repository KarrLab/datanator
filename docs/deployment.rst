Deployment
============
The following instructions describe how to deploy ``datanator`` to the heroku server

We are deploying the backend API server via a container using the karrlab/wc_env_dependencies:latest image.

The commands for deploying the container are the following::

  heroku login
  heroku container:login
  heroku container:push web -a datanator
  heroku container:release web -a datanator

In order to change the configuration of the container, look at the Dockerfile for datanator. The gunicorn production server can be
adjusted accordingly in order to accommodate the number of users. 


Contact `Saahith <mailto:saahith116@gmail.com>`_ for any questions regarding deployment
