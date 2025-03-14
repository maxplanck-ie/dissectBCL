import yaml
from urllib.parse import urlparse
import requests

class Test_Contam():
    def test_links(self):
        # contam is a function that is not yet implemented
        with open("contaminome.yml") as f:
            genomes = yaml.safe_load(f)
        _links = []
        for gtype in genomes:
            for genome in genomes[gtype]:
                _links.append(
                    genomes[gtype][genome]['URL']
                )
        for _link in _links:
            _lp = urlparse(_link)
            # https link
            assert _lp.scheme == 'https', f"URL is not an https link: {_link}"
            # Domain not empty
            assert _lp.netloc, f"URL seems to not have a domain: {_link}"
            # responsive
            response = requests.head(_link, allow_redirects=True, timeout=10)
            assert response.status_code == 200, f"URL not reachable (status {response.status_code}): {_link}"
    
    def test_failure(self):
        assert False