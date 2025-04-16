# Configuration

BioV can be configured through:

1. Environment variables
2. `.env` files
3. Runtime settings

## Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `BIOV_HOME` | Platform cache dir | Custom cache directory |
| `BIOV_CACHE_HTTP` | `True` | Enable HTTP caching |

```python
from biov.config import settings
print(settings.cache_http)
```
