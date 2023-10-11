/**
 * Helper functions for the embedded graph colors
 */
import * as d3 from "d3";
import { interpolateRainbow, interpolateCool } from "d3-scale-chromatic";
import memoize from "memoize-one";
import * as globals from "../../globals";
import parseRGB from "../parseRGB";
import { range } from "../range";
import {
  Dataframe,
  DataframeColumn,
  DataframeValueArray,
  LabelType,
} from "../dataframe";
import {
  Field,
  Schema,
  CategoricalAnnotationColumnSchema,
} from "../../common/types/schema";
import { Query } from "../../annoMatrix/query";
import { Genesets } from "../../reducers/genesets";
import {
  ColorsState,
  ConvertedUserColors,
  UserColor,
} from "../../reducers/colors";
import { isDataframeDictEncodedColumn } from "../dataframe/types";

/**
 * Given a color mode & accessor, generate an annoMatrix query that will
 * fulfill it
 */
export function createColorQuery(
  colorMode: ColorsState["colorMode"],
  colorByAccessor: string,
  schema: Schema,
  genesets: Genesets,
): [Field, Query] | null {
  if (!colorMode || !colorByAccessor || !schema || !genesets) return null;

  switch (colorMode) {
    case "color by categorical metadata":
    case "color by continuous metadata": {
      return [Field.obs, colorByAccessor];
    }
    case "color by expression": {
      const varIndex = schema?.annotations?.var?.index;
      if (!varIndex) return null;
      return [
        Field.X,
        {
          where: {
            field: "var",
            column: varIndex,
            value: colorByAccessor,
          },
        },
      ];
    }
    case "color by geneset mean expression": {
      const varIndex = schema?.annotations?.var?.index;

      if (!varIndex) return null;
      if (!genesets) return null;

      const geneset = genesets.get(colorByAccessor);
      const setGenes = [...(geneset?.genes.keys() || [])];

      return [
        Field.X,
        {
          summarize: {
            method: "mean",
            field: "var",
            column: varIndex,
            values: setGenes,
          },
        },
      ];
    }
    default: {
      return null;
    }
  }
}

function _defaultColors(nObs: number): ColorTable {
  return {
    rgb: new Float32Array(3 * nObs).fill(0),
    scale: undefined,
  };
}

const defaultColors = memoize(_defaultColors);

export type ColorRange =
  | d3.ScaleSequentialBase<string, never>["range"]
  | d3.ScaleQuantile<number, never>["range"];
export interface ColorTable {
  rgb: Float32Array;
  scale:
    | {
        (idx: number): d3.RGBColor;
        range?(): ColorRange;
        domain(): number[];
      }
    | d3.ScaleQuantile<number, never>
    | d3.ScaleSequential<string, never>
    | undefined;
}
/**
 * Create colors scale and RGB array and return as object.
 *
 * @param colorMode - categorical, etc.
 * @param colorByAccessor - the annotation label name
 * @param colorByDataframe - the actual color-by data
 * @param schema - the entire schema
 * @param userColors - optional user color table
 * @returns colors scale and RGB array as an object
 */
function _createColorTable(
  colorMode: string | null,
  colorByAccessor: LabelType | null,
  colorByData: Dataframe | null,
  schema: Schema,
  userColors: ConvertedUserColors | null = null,
): ColorTable {
  if (colorMode === null || colorByData === null) {
    return defaultColors(schema.dataframe.nObs);
  }

  switch (colorMode) {
    case "color by categorical metadata": {
      if (colorByAccessor === null) return defaultColors(schema.dataframe.nObs);

      const col = colorByData.col(colorByAccessor);

      if (userColors && colorByAccessor in userColors) {
        return createUserColors(col, colorByAccessor, schema, userColors);
      }

      return createColorsByCategoricalMetadata(col, colorByAccessor, schema);
    }
    case "color by continuous metadata": {
      if (colorByAccessor === null) return defaultColors(schema.dataframe.nObs);
      const col = colorByData.col(colorByAccessor);
      const { min, max } = col.summarizeContinuous();
      return createColorsByContinuousMetadata(col.asArray(), min, max);
    }
    case "color by expression": {
      const col = colorByData.icol(0);
      const { min, max } = col.summarizeContinuous();
      return createColorsByContinuousMetadata(col.asArray(), min, max);
    }
    case "color by geneset mean expression": {
      const col = colorByData.icol(0);
      const { min, max } = col.summarizeContinuous();
      return createColorsByContinuousMetadata(col.asArray(), min, max);
    }
    default: {
      return defaultColors(schema.dataframe.nObs);
    }
  }
}
export const createColorTable = memoize(_createColorTable);

/**
 * Create two category label-indexed objects:
 *  - colors: maps label to RGB triplet for that label (used by graph, etc)
 *  - scale: function which given label returns d3 color scale for label
 * Order doesn't matter - everything is keyed by label value.
 */
export function loadUserColorConfig(userColors: {
  [category: string]: { [label: string]: string };
}): ConvertedUserColors {
  const convertedUserColors = {} as ConvertedUserColors;

  Object.keys(userColors).forEach((category) => {
    const [colors, scaleMap] = Object.keys(userColors[category]).reduce(
      (acc, label) => {
        const color = parseRGB(userColors[category][label]);

        acc[0][label] = color;

        acc[1][label] = d3.rgb(255 * color[0], 255 * color[1], 255 * color[2]);
        return acc;
      },
      [{} as UserColor["colors"], {} as { [label: string]: d3.RGBColor }],
    );

    const scale = (label: string) => scaleMap[label];

    const userColor = { colors, scale };

    convertedUserColors[category] = userColor;
  });

  return convertedUserColors;
}

function _createUserColors(
  col: DataframeColumn,
  colorAccessor: LabelType,
  schema: Schema,
  userColors: ConvertedUserColors,
): ColorTable {
  const data = col.asArray();
  const { colors, scale: scaleByLabel } = userColors[colorAccessor];
  let newColors = colors;
  // convert the keys in the color dictionary defined by the schema
  // to their corresponding codes.
  if (isDataframeDictEncodedColumn(col)) {
    newColors = Object.fromEntries(
      Object.entries(colors).map((row) => [col.invCodeMapping[row[0]], row[1]]),
    );
  }
  const rgb = createRgbArray(data, newColors);

  // color scale function param is INDEX (offset) into schema categories. It is NOT label value.
  // See createColorsByCategoricalMetadata() for another example.
  // TODO: #35 Use type guards to insure type instead of casting
  let { categories } = schema.annotations.obsByName[
    colorAccessor
  ] as CategoricalAnnotationColumnSchema;
  if (isDataframeDictEncodedColumn(col)) {
    categories = categories.map((cat) => col.codeMapping[cat as number]);
  }
  const categoryMap = new Map();

  categories?.forEach((label, idx) => categoryMap.set(idx, label));
  const scale = (idx: number) => scaleByLabel(categoryMap.get(idx));
  scale.domain = () => [0, 0];
  return { rgb, scale };
}
const createUserColors = memoize(_createUserColors);

interface CategoryColors {
  [category: string]: [number, number, number];
}

function _createColorsByCategoricalMetadata(
  col: DataframeColumn,
  colorAccessor: LabelType,
  schema: Schema,
): ColorTable {
  const data = col.asArray();
  // TODO: #35 Use type guards to insure type instead of casting
  const { categories } = schema.annotations.obsByName[
    colorAccessor
  ] as CategoricalAnnotationColumnSchema;

  const scale = d3
    .scaleSequential(interpolateRainbow)
    .domain([0, categories?.length || 0]);

  /* pre-create colors - much faster than doing it for each obs */
  const colors = categories?.reduce((acc: CategoryColors, cat, idx) => {
    acc[cat as string] = parseRGB(scale(idx));
    return acc;
  }, {});
  const rgb = createRgbArray(data, colors);
  return { rgb, scale };
}

const createColorsByCategoricalMetadata = memoize(
  _createColorsByCategoricalMetadata,
);

function createRgbArray(data: DataframeValueArray, colors?: CategoryColors) {
  const rgb = new Float32Array(3 * data.length);

  if (!colors) {
    throw new Error("`colors` is undefined");
  }

  for (let i = 0, len = data.length; i < len; i += 1) {
    const label = data[i];
    rgb.set(colors[String(label)], 3 * i);
  }

  return rgb;
}

function _createColorsByContinuousMetadata(
  data: DataframeValueArray,
  min: number,
  max: number,
): ColorTable {
  const colorBins = 100;
  const scale = d3
    .scaleQuantile()
    .domain([min, max])
    .range(range(colorBins - 1, -1, -1));

  /* pre-create colors - much faster than doing it for each obs */
  const colors = new Array(colorBins);
  for (let i = 0; i < colorBins; i += 1) {
    colors[i] = parseRGB(interpolateCool(i / colorBins));
  }

  const nonFiniteColor = parseRGB(globals.nonFiniteCellColor);

  const rgb = new Float32Array(3 * data.length);
  for (let i = 0, len = data.length; i < len; i += 1) {
    const val = data[i];

    if (Number.isFinite(val)) {
      const c = scale(val as number);
      rgb.set(colors[c], 3 * i);
    } else {
      rgb.set(nonFiniteColor, 3 * i);
    }
  }
  return { rgb, scale };
}

export const createColorsByContinuousMetadata = memoize(
  _createColorsByContinuousMetadata,
);
