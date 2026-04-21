from cache.directory import CacheDirectory, CacheItem, get_cache_item

test_item = CacheItem(
    name="Test Item",
    cached="https://drive.google.com/uc?id=",
    unpinned="",
    pinned=""
)

custom_directory: CacheDirectory = {
    "test": {
        "alias": "v1",
        "v1": test_item
    }
}

def test_directory_fetch():
    assert get_cache_item(("test", "v1"), custom_directory=custom_directory) == test_item
    assert get_cache_item(("test", "alias"), custom_directory=custom_directory) == test_item
