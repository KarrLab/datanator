FROM node:8-alpine

WORKDIR /app

RUN npm i npm@latest -g

COPY . /app

EXPOSE 3001
