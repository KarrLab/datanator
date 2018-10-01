from flask import Flask, render_template

import json



def render(schema):
    app = Flask(__name__)
    @app.route('/')
    def hello_world(variable = ""):
        print(schema)
        return render_template('index2.html', variable=schema)
    app.run().hello_world(schema)


if __name__ == '__main__':
    schema = {
"title": "Person",
"type": "object",
"properties": {
"name": {
  "type": "string",
  "description": "First and Last name",
  "minLength": 4,
  "default": "Jeremy Dorn"
}}}
    render(schema)

    