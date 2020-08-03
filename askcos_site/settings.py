"""
Django settings for askcos_site project.

For more information on this file, see
https://docs.djangoproject.com/en/2.2/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/2.2/ref/settings/
"""

import os

# Get settings from separate celeryconfig.py
from askcos_site.askcos_celery.celeryconfig import *


PROJECT_PATH = os.path.dirname(__file__)
BASE_DIR = os.path.dirname(PROJECT_PATH)

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = 'notsosecret'
# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = os.getenv('DJANGO_DEBUG', False)

ALLOWED_HOSTS = ['0.0.0.0']
if os.environ.get('CURRENT_HOST'):
    ALLOWED_HOSTS.append(os.environ.get('CURRENT_HOST'))

TEMPLATE_LOADERS = [
    'django.template.loaders.filesystem.Loader',
    'django.template.loaders.app_directories.Loader'
]

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [
            os.path.join(PROJECT_PATH, 'templates'),
        ],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                # Insert your TEMPLATE_CONTEXT_PROCESSORS here or use this
                # list if you haven't customized them:
                'askcos_site.processors.customization',
                'django.contrib.auth.context_processors.auth',
                'django.template.context_processors.debug',
                'django.template.context_processors.i18n',
                'django.template.context_processors.media',
                'django.template.context_processors.static',
                'django.template.context_processors.tz',
                'django.template.context_processors.request',
                'django.contrib.messages.context_processors.messages',
            ],
        },
    },
]

INSTALLED_APPS = (
    'askcos_site.main',
    'django.contrib.auth',
    'django.contrib.admin',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'django_extensions',
    'rest_framework',
    'rest_framework_jwt',
)

MIDDLEWARE = (
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
)

ROOT_URLCONF = 'askcos_site.urls'
WSGI_APPLICATION = 'askcos_site.wsgi.application'

# Base django auth settings
LOGIN_URL = '/registration/login'
LOGIN_REDIRECT_URL = '/'

# django-registration settings
ACCOUNT_ACTIVATION_DAYS = 7
REGISTRATION_OPEN = True
REGISTRATION_SALT = 'saltystring'  # does not need to be secret, for namespace

# Email address for automated correspondence
DEFAULT_FROM_EMAIL = 'no-reply@askcos.mit.edu'

# Where are user settings / banlists / etc. saved?
# NOTE: we recommend relocating the db to an ssd for speed
DATABASES = {'default': {
    'ENGINE': 'django.db.backends.mysql',
    'NAME': os.getenv('MYSQL_DATABASE', 'askcos_db'),
    'USER': os.getenv('MYSQL_USER', 'root'),
    'PASSWORD': os.getenv('MYSQL_ROOT_PASSWORD', 'password'),
    'HOST': os.getenv('MYSQL_HOST', 'mysql'),
    'PORT': '3306',
}}

# django rest framework settings
REST_FRAMEWORK = {
    'DEFAULT_AUTHENTICATION_CLASSES': [
        'rest_framework.authentication.BasicAuthentication',
        'rest_framework.authentication.SessionAuthentication',
        'rest_framework_jwt.authentication.JSONWebTokenAuthentication',
    ]
}

LANGUAGE_CODE = 'en-us'
TIME_ZONE = 'UTC'
USE_I18N = True
USE_L10N = True
USE_TZ = True

# Static files (CSS, JavaScript, Images)
STATIC_ROOT = os.path.join(os.path.dirname(PROJECT_PATH), 'static')  # path that collectstatic places static files
STATIC_URL = os.getenv('STATIC_URL', '/static/')  # url from which static files are served
STATICFILES_DIRS = (
    os.path.join(PROJECT_PATH, 'static'),
)  # directories containing static files to be moved copied to STATIC_ROOT during deployment

# Media files
MEDIA_ROOT = os.path.join(os.path.dirname(PROJECT_PATH), 'media')  # location for placing user uploaded files
MEDIA_URL = '/media/'  # url from which user uploaded files are served

################################################################################
# Define databases to replicate Make-It settings
################################################################################

LOCAL_STORAGE = {}


# Try importing any custom settings here so they replace the defaults above
try:
    from askcos_site.custom_settings import *
except ImportError:
    pass
