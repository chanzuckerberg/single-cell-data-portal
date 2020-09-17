type Replace = { [target: string]: string };

export function apiTemplateToUrl(template: string, replace: Replace) {
  let result = template;

  for (const [target, value] of Object.entries(replace)) {
    if (result.indexOf(target) < 0) {
      throw Error("Cannot find target");
    }

    result = result.replace(`{${target}}`, value);
  }

  return result;
}
