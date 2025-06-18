import { extractDatasetType } from "./outputStyle";

export function getDataFiles() {
    const dataFiles = import.meta.glob('../../public/data/output/**', { query: '?raw' });
    return Object.keys(dataFiles)
        .map(path => path.substring("../../public/data/output/".length))
}

export function getDatasets() {
    const files = getDataFiles();
    return files.filter(file => file.startsWith("logs/datasets-"))
        .map(file => file.substring("logs/datasets-".length))
        .map(file => file.slice(0, -".yaml".length))
        .map(file => extractDatasetType(file))
}
