import { parse } from 'yaml'

import dmmmYaml from '../../public/data/configs/dmmm.yaml?raw'
import praYaml from '../../public/data/configs/dmmm.yaml?raw'

const configs: Record<string, Record<string, unknown>> = {
    'dmmm': parse(dmmmYaml),
    'pra': parse(praYaml),
}

const datasets = Object.entries(configs).map((([type, entry]) =>
    (entry["datasets"] as Record<string, unknown>[]).map(dataset => ({ ...dataset, type }))
)).flat()

function getDatasetData(datasetName: string): Dataset {

}