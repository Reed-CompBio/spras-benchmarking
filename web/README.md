# web

This module is an [Astro](https://astro.build/) project which wraps the output from SPRAS
into a presentable webpage. See the output: https://reed-compbio.github.io/spras-benchmarking/

## Building

To build this, you need [`pnpm`](https://pnpm.io/). It is recommended to use a node version manager
([nvm](https://github.com/nvm-sh/nvm) for mac/linux, [nvm-windows](https://github.com/coreybutler/nvm-windows) for windows),
to install `nodejs` and `npm` (at the time of writing, this would be node `v22`), and use `npm` to install `pnpm`:

```sh
npm install --global pnpm
```

After this, you can install the dependencies (make sure your current working directory is `web`):

```sh
pnpm install
```

Then, assuming your data is in `public/data`, build the website:

```sh
pnpm run build
```
