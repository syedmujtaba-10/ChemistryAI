declare module "rdkit" {
    interface RDKitModule {
      get_mol(smiles: string): any;
      get_svg(): string;
      draw_chemical_structure(): string;
    }
  
    function initRDKitModule(): Promise<RDKitModule>;
  
    export { initRDKitModule, RDKitModule };
  }
  