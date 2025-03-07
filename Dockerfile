# Usamos una imagen base de Python 3.10
FROM python:3.10-slim-buster

# Establecemos el directorio de trabajo en el contenedor
WORKDIR /usr/src/app

# Copiamos el archivo requirements.txt al directorio de trabajo
COPY requirements.txt ./

# Instalamos las dependencias
RUN pip install --no-cache-dir -r requirements.txt

# Copiamos el resto del código de la aplicación al directorio de trabajo
COPY app/ ./app/

# Exponemos el puerto en el que se ejecutará la aplicación
EXPOSE 8080

# Especificamos el comando para iniciar la aplicación
CMD ["python", "./montu-app/app.py"]
