from google_auth_oauthlib import flow
import google.auth

def login():
    '''Run the following to generate a default credentials file
    
    $ gcloud auth application-default login
    
    https://google-auth.readthedocs.io/en/latest/reference/google.auth.html#google.auth.default.'''
    
    credentials, gcp_project = google.auth.default()
    return credentials, gcp_project


def active_login(client_secrets):
    '''# TODO: Uncomment the line below to set the `launch_browser` variable.
    # launch_browser = True
    #
    # The `launch_browser` boolean variable indicates if a local server is used
    # as the callback URL in the auth flow. A value of `True` is recommended,
    # but a local server does not work if accessing the application remotely,
    # such as over SSH or from a remote Jupyter notebook.'''
    appflow = flow.InstalledAppFlow.from_client_secrets_file(
        client_secrets, scopes=["https://www.googleapis.com/auth/bigquery"]
    )
    launch_browser = False
    if launch_browser:
        appflow.run_local_server()
    else:
        appflow.run_console()
    credentials = appflow.credentials
    return credentials