FROM python:2.7

RUN pip install \
    'numpy==1.16.6' \
    'scikit-learn==0.17.1' \
    'scipy==1.2.3' \
    'matplotlib==2.2.5' \
    'biopython==1.76' \
    'azimuth==2.0'

ENTRYPOINT ["python"]
