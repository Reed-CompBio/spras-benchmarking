export function getPaths() {
    const dataFiles = import.meta.glob('../../public/data/output/**', { query: '?raw' });
    return Object.keys(dataFiles)
        .map(path => path.substring("../../public/data/output/".length))
}