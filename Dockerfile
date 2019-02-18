# Use an official Python runtime as a parent image
FROM python:3.7-slim

# Set the working directory to /app
WORKDIR /app

# Copy the current directory contents into the container at /app
COPY *.py requirements.txt /app/

# Install any needed packages specified in requirements.txt
RUN pip install --trusted-host pypi.python.org -r requirements.txt

RUN apt-get update && apt-get install -y tk

# Make port 80 available to the world outside this container
# EXPOSE 80

# Define environment variable
# ENV NAME World
ENV DISPLAY 127.0.0.1:0

# Run app.py when the container launches
CMD ["python", "one_tissue_analysis.py"]
