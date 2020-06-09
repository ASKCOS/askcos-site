import os


def customization(request):
    customizable = {
        'organization': os.environ.get('ORGANIZATION', 'MIT'),
        'contact_email': os.environ.get('CONTACT_EMAIL', 'mlpds_support@mit.edu'),
        'support_emails': os.environ.get('SUPPORT_EMAILS', 'incoming+mlpds-mit-askcos-askcos-12564933-issue-@incoming.gitlab.com,mlpds_support@mit.edu'),
        'email_enabled': bool(os.environ.get('EMAIL_ENABLED', False)),
        'version_number': os.environ.get('VERSION_NUMBER', '0.x'),
        'update_date': os.environ.get('UPDATE_DATE', ''),
    }
    return customizable
