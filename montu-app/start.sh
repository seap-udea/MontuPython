#!/bin/bash
gunicorn -w 2 app:server — bind=0.0.0.0:8050
