FROM mongo

# Create app directory
WORKDIR /usr/src/configs

# Install app dependencies
COPY mongoSetup.js .
COPY setup.sh .

RUN ["chmod", "+x", "./setup.sh"]
CMD ["./setup.sh"]