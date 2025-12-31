interface Output {
    dataType: string;
    datasetCategory: string;
    datasetName: string;
    algorithm: string;
    paramsHash: string;
}

function extractPrefix(name: string, prefixName: string, prefixes: string[]): { prefix: string, name: string } {
    const foundPrefix = prefixes.find(prefix => name.startsWith(prefix))
    
    if (!foundPrefix) {
        throw new Error(`${name} should begin with a ${prefixName} (one of ${prefixes})`)
    }

    return { 
        prefix: foundPrefix,
        name: name.substring(foundPrefix.length)
    }
}

const dataTypes = [
    'pra',
    'dmmm'
]

export function extractDatasetType(name: string): { type: string, name: string } {
    const { prefix, name: newName } = extractPrefix(name, "dataset type", dataTypes)
    return { type: prefix, name: newName }
}

const dataCategories = {
    'diseases': {
        name: 'DISEASES',
        directory: 'diseases'
    },
    'depmap': {
        name: 'DepMap',
        directory: 'depmap'
    },
    'hiv': {
        name: 'HIV',
        directory: 'hiv'
    },
    'rn': {
        name: 'ResponseNet',
        directory: 'rn-muscle-skeletal'
    },
    'yeast': {
        name: 'Yeast',
        directory: 'yeast-osmotic-stress'
    },
}

// TODO: replace this once we have proper dataset categories
export function extractDatasetCategory(name: string): { category: string, name: string } {
    const { prefix, name: newName } = extractPrefix(name, "dataset category", Object.keys(dataCategories))
    return { category: prefix, name: newName.slice(1) }
}

export function parseOutputString(str: string): Output {
    const components = str.split("-");
    let dataType;
    let datasetCategory;
    let datasetName;
    let algorithm;
    let paramsHash;

    if (components.length === 5) {
        // This is a slug URL (type-...)
        [dataType, datasetCategory, datasetName, algorithm, paramsHash] = components
    } else if (components.length === 3) {
        // This is fetched straight from the folder - we ignored -params previously
        [datasetName, algorithm, paramsHash] = components
    } else {
        throw new Error(`Unexpected length of components in ${components}.`)
    }


    // We didn't get a data type in the first passthrough - lets extract the data
    // type from the name
    if (!dataType || !datasetCategory) {
        const { type, name: name1 } = extractDatasetType(datasetName);
        const { category, name } = extractDatasetCategory(name1)
        dataType = type;
        datasetCategory = category;
        datasetName = name;
    }

    return {
        dataType,
        datasetCategory,
        datasetName,
        algorithm,
        paramsHash
    }
}

export function styleOutput(output: Output): string {
    return `${output.dataType}-${output.datasetCategory}-${output.datasetName}-${output.algorithm}-${output.paramsHash}`
}

export function asFolderName(output: Output): string {
    return `${output.dataType}${output.datasetCategory}_${output.datasetName}-${output.algorithm}-params-${output.paramsHash}`
}

export function algorithmDocumentationUrl(algorithm: string): string {
    const map: Record<string, string> = {
        "omicsintegrator1": "oi1",
        "omicsintegrator2": "oi2"
    }

    const foundAlgorithm = algorithm in map ? map[algorithm] : algorithm
    return `https://spras.readthedocs.io/en/latest/prms/${foundAlgorithm}.html`;
}
